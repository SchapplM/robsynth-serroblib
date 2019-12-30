% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
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
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(7);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t15 = qJ(1) + pkin(7);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:06
	% EndTime: 2019-12-29 15:46:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (41->15), mult. (52->20), div. (0->0), fcn. (36->6), ass. (0->9)
	t102 = r_i_i_C(2) + qJ(3);
	t97 = sin(pkin(8));
	t101 = qJD(1) * t97;
	t100 = qJD(4) * t97;
	t99 = -pkin(2) + (-pkin(3) - r_i_i_C(1)) * cos(pkin(8)) + (-r_i_i_C(3) - qJ(4)) * t97;
	t96 = qJ(1) + pkin(7);
	t95 = cos(t96);
	t94 = sin(t96);
	t1 = [-t94 * t100 + t95 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t102 * t94 + t99 * t95) * qJD(1), 0, qJD(1) * t95, -t94 * t101, 0; t95 * t100 + t94 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t102 * t95 + t99 * t94) * qJD(1), 0, qJD(1) * t94, t95 * t101, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (106->27), mult. (164->47), div. (0->0), fcn. (138->8), ass. (0->22)
	t40 = sin(pkin(8));
	t53 = qJD(1) * t40;
	t52 = qJD(4) * t40;
	t51 = pkin(6) + r_i_i_C(3) - qJ(3);
	t41 = cos(pkin(8));
	t42 = sin(qJ(5));
	t43 = cos(qJ(5));
	t50 = -t40 * t43 + t41 * t42;
	t49 = t40 * t42 + t41 * t43;
	t39 = qJ(1) + pkin(7);
	t38 = cos(t39);
	t48 = t49 * t38;
	t37 = sin(t39);
	t47 = t49 * t37;
	t46 = qJD(1) * t50;
	t45 = t50 * qJD(5);
	t44 = -qJ(4) * t40 - pkin(2) + (-pkin(3) - pkin(4)) * t41;
	t36 = -qJD(1) * t48 + t37 * t45;
	t35 = qJD(5) * t47 + t38 * t46;
	t34 = qJD(1) * t47 + t38 * t45;
	t33 = -qJD(5) * t48 + t37 * t46;
	t1 = [-t37 * t52 + t36 * r_i_i_C(1) + t35 * r_i_i_C(2) + t38 * qJD(3) + (-cos(qJ(1)) * pkin(1) + t51 * t37 + t44 * t38) * qJD(1), 0, qJD(1) * t38, -t37 * t53, t33 * r_i_i_C(1) + t34 * r_i_i_C(2); t38 * t52 - t34 * r_i_i_C(1) + t33 * r_i_i_C(2) + t37 * qJD(3) + (-sin(qJ(1)) * pkin(1) - t51 * t38 + t44 * t37) * qJD(1), 0, qJD(1) * t37, t38 * t53, -t35 * r_i_i_C(1) + t36 * r_i_i_C(2); 0, 0, 0, 0, (t50 * r_i_i_C(1) + t49 * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end