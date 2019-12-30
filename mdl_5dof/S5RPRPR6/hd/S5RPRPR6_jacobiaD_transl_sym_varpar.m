% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:27
	% EndTime: 2019-12-29 16:46:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:22
	% EndTime: 2019-12-29 16:46:22
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:27
	% EndTime: 2019-12-29 16:46:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(8);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:22
	% EndTime: 2019-12-29 16:46:22
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t42 = qJ(1) + pkin(8);
	t40 = qJ(3) + t42;
	t38 = sin(t40);
	t39 = cos(t40);
	t41 = qJD(1) + qJD(3);
	t44 = (-r_i_i_C(1) * t39 + r_i_i_C(2) * t38) * t41;
	t43 = (-r_i_i_C(1) * t38 - r_i_i_C(2) * t39) * t41;
	t1 = [(-cos(t42) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t44, 0, t44, 0, 0; t43 + (-sin(t42) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t43, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:22
	% EndTime: 2019-12-29 16:46:22
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (84->14), mult. (46->16), div. (0->0), fcn. (26->6), ass. (0->11)
	t33 = r_i_i_C(3) + qJ(4);
	t28 = qJ(1) + pkin(8);
	t26 = qJ(3) + t28;
	t24 = sin(t26);
	t27 = qJD(1) + qJD(3);
	t32 = t27 * t24;
	t25 = cos(t26);
	t31 = t27 * t25;
	t30 = t24 * qJD(4) + (-pkin(3) + r_i_i_C(2)) * t32 + t33 * t31;
	t29 = r_i_i_C(2) * t31 + t25 * qJD(4) + (-pkin(3) * t25 - t33 * t24) * t27;
	t1 = [(-cos(t28) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t29, 0, t29, t31, 0; (-sin(t28) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t30, 0, t30, t32, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:46:27
	% EndTime: 2019-12-29 16:46:27
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (167->26), mult. (122->41), div. (0->0), fcn. (76->8), ass. (0->20)
	t47 = qJ(1) + pkin(8);
	t45 = qJ(3) + t47;
	t43 = sin(t45);
	t44 = cos(t45);
	t49 = cos(qJ(5));
	t58 = qJD(5) * t49;
	t46 = qJD(1) + qJD(3);
	t48 = sin(qJ(5));
	t61 = t46 * t48;
	t63 = t43 * t58 + t44 * t61;
	t62 = t46 * t44;
	t60 = t46 * t49;
	t59 = qJD(5) * t48;
	t57 = -pkin(3) - pkin(7) - r_i_i_C(3);
	t55 = t44 * t60;
	t53 = t44 * t59;
	t52 = t44 * t58;
	t51 = r_i_i_C(2) * t55 + qJ(4) * t62 + (-r_i_i_C(2) * t59 + t57 * t46 + qJD(4)) * t43 + t63 * r_i_i_C(1);
	t50 = -r_i_i_C(2) * t53 + r_i_i_C(1) * t52 + t44 * qJD(4) + (t57 * t44 + (-r_i_i_C(1) * t48 - r_i_i_C(2) * t49 - qJ(4)) * t43) * t46;
	t1 = [(-cos(t47) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t50, 0, t50, t62, -t63 * r_i_i_C(2) + (-t43 * t59 + t55) * r_i_i_C(1); (-sin(t47) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t51, 0, t51, t46 * t43, (-t43 * t61 + t52) * r_i_i_C(2) + (t43 * t60 + t53) * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t49 + r_i_i_C(2) * t48) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end