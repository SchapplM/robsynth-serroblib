% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t19 = pkin(1) * qJD(1);
	t16 = qJ(1) + qJ(2);
	t13 = sin(t16);
	t14 = cos(t16);
	t15 = qJD(1) + qJD(2);
	t18 = (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t15;
	t17 = (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t15;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t19 + t17, t17, 0, 0, 0; cos(qJ(1)) * t19 + t18, t18, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t25 = pkin(1) * qJD(1);
	t22 = qJ(1) + qJ(2);
	t19 = pkin(8) + t22;
	t17 = sin(t19);
	t18 = cos(t19);
	t21 = qJD(1) + qJD(2);
	t24 = (cos(t22) * pkin(2) + r_i_i_C(1) * t18 - r_i_i_C(2) * t17) * t21;
	t23 = (-pkin(2) * sin(t22) - r_i_i_C(1) * t17 - r_i_i_C(2) * t18) * t21;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t25 + t23, t23, 0, 0, 0; cos(qJ(1)) * t25 + t24, t24, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (109->16), mult. (66->20), div. (0->0), fcn. (40->8), ass. (0->15)
	t69 = qJD(1) + qJD(2);
	t83 = pkin(2) * t69;
	t82 = pkin(3) + r_i_i_C(1) * cos(pkin(9));
	t81 = -qJD(4) - r_i_i_C(2) * t69 * sin(pkin(9));
	t70 = qJ(1) + qJ(2);
	t67 = pkin(8) + t70;
	t65 = sin(t67);
	t79 = t69 * t65;
	t66 = cos(t67);
	t78 = t69 * t66;
	t77 = pkin(1) * qJD(1);
	t76 = qJ(4) * t69;
	t74 = t65 * t76 + r_i_i_C(3) * t79 + cos(t70) * t83 + t81 * t66 + t82 * t78;
	t73 = -sin(t70) * t83 + r_i_i_C(3) * t78 + t66 * t76 + (-t69 * t82 - t81) * t65;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t77 + t73, t73, 0, t79, 0; cos(qJ(1)) * t77 + t74, t74, 0, -t78, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (194->30), mult. (122->43), div. (0->0), fcn. (78->9), ass. (0->23)
	t117 = pkin(1) * qJD(1);
	t103 = qJD(1) + qJD(2);
	t104 = qJ(1) + qJ(2);
	t100 = pkin(8) + t104;
	t95 = sin(t100);
	t116 = t103 * t95;
	t96 = cos(t100);
	t115 = t103 * t96;
	t102 = pkin(9) + qJ(5);
	t98 = sin(t102);
	t114 = t103 * t98;
	t99 = cos(t102);
	t113 = t103 * t99;
	t112 = qJD(5) * t98;
	t111 = qJD(5) * t99;
	t110 = t95 * t114;
	t109 = t96 * t113;
	t108 = (-r_i_i_C(1) * t98 - r_i_i_C(2) * t99) * qJD(5);
	t105 = -pkin(7) - qJ(4);
	t97 = cos(pkin(9)) * pkin(4) + pkin(3);
	t107 = t97 * t115 + r_i_i_C(1) * t109 + r_i_i_C(3) * t116 + pkin(2) * t103 * cos(t104) + (-r_i_i_C(2) * t114 - qJD(4)) * t96 + (-t103 * t105 + t108) * t95;
	t106 = r_i_i_C(2) * t110 + r_i_i_C(3) * t115 + t95 * qJD(4) + t96 * t108 + (-pkin(2) * sin(t104) - t105 * t96 + (-r_i_i_C(1) * t99 - t97) * t95) * t103;
	t1 = [0, 0, 0, 0, t108; -sin(qJ(1)) * t117 + t106, t106, 0, t116, (t95 * t112 - t109) * r_i_i_C(2) + (-t95 * t111 - t96 * t114) * r_i_i_C(1); cos(qJ(1)) * t117 + t107, t107, 0, -t115, (-t96 * t112 - t95 * t113) * r_i_i_C(2) + (t96 * t111 - t110) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end