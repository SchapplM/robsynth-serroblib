% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:11:40
	% EndTime: 2019-12-29 18:11:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:11:35
	% EndTime: 2019-12-29 18:11:35
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
	% StartTime: 2019-12-29 18:11:35
	% EndTime: 2019-12-29 18:11:35
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:11:35
	% EndTime: 2019-12-29 18:11:35
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t49 = pkin(1) * qJD(1);
	t46 = qJ(1) + qJ(2);
	t42 = pkin(8) + t46;
	t40 = sin(t42);
	t41 = cos(t42);
	t45 = qJD(1) + qJD(2);
	t48 = (-pkin(2) * cos(t46) - r_i_i_C(1) * t41 + t40 * r_i_i_C(2)) * t45;
	t47 = (-pkin(2) * sin(t46) - r_i_i_C(1) * t40 - r_i_i_C(2) * t41) * t45;
	t1 = [-cos(qJ(1)) * t49 + t48, t48, 0, 0, 0; -sin(qJ(1)) * t49 + t47, t47, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:11:40
	% EndTime: 2019-12-29 18:11:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (92->14), mult. (50->17), div. (0->0), fcn. (28->6), ass. (0->12)
	t38 = r_i_i_C(3) + qJ(4);
	t32 = qJ(1) + qJ(2);
	t28 = pkin(8) + t32;
	t26 = sin(t28);
	t31 = qJD(1) + qJD(2);
	t37 = t31 * t26;
	t27 = cos(t28);
	t36 = t31 * t27;
	t35 = pkin(1) * qJD(1);
	t34 = r_i_i_C(2) * t37 + t26 * qJD(4) + (-pkin(2) * sin(t32) - pkin(3) * t26) * t31 + t38 * t36;
	t33 = r_i_i_C(2) * t36 + t27 * qJD(4) + (-pkin(2) * cos(t32) - pkin(3) * t27 - t38 * t26) * t31;
	t1 = [-cos(qJ(1)) * t35 + t33, t33, 0, t36, 0; -sin(qJ(1)) * t35 + t34, t34, 0, t37, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:11:40
	% EndTime: 2019-12-29 18:11:40
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (175->26), mult. (126->41), div. (0->0), fcn. (78->8), ass. (0->21)
	t51 = qJ(1) + qJ(2);
	t47 = pkin(8) + t51;
	t45 = sin(t47);
	t46 = cos(t47);
	t53 = cos(qJ(5));
	t62 = qJD(5) * t53;
	t50 = qJD(1) + qJD(2);
	t52 = sin(qJ(5));
	t66 = t50 * t52;
	t68 = t45 * t62 + t46 * t66;
	t67 = t50 * t46;
	t65 = t50 * t53;
	t64 = pkin(1) * qJD(1);
	t63 = qJD(5) * t52;
	t61 = -pkin(3) - pkin(7) - r_i_i_C(3);
	t59 = t46 * t65;
	t57 = t46 * t63;
	t56 = t46 * t62;
	t55 = -pkin(2) * t50 * sin(t51) + r_i_i_C(2) * t59 + qJ(4) * t67 + (-r_i_i_C(2) * t63 + t61 * t50 + qJD(4)) * t45 + t68 * r_i_i_C(1);
	t54 = -r_i_i_C(2) * t57 + r_i_i_C(1) * t56 + t46 * qJD(4) + (-pkin(2) * cos(t51) + t61 * t46 + (-r_i_i_C(1) * t52 - r_i_i_C(2) * t53 - qJ(4)) * t45) * t50;
	t1 = [-cos(qJ(1)) * t64 + t54, t54, 0, t67, -t68 * r_i_i_C(2) + (-t45 * t63 + t59) * r_i_i_C(1); -sin(qJ(1)) * t64 + t55, t55, 0, t50 * t45, (-t45 * t66 + t56) * r_i_i_C(2) + (t45 * t65 + t57) * r_i_i_C(1); 0, 0, 0, 0, (-r_i_i_C(1) * t53 + r_i_i_C(2) * t52) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end