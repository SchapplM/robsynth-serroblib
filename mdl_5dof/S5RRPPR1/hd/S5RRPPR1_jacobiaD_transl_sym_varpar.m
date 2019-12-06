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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 18:19:03
	% EndTime: 2019-12-05 18:19:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:19:03
	% EndTime: 2019-12-05 18:19:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:19:03
	% EndTime: 2019-12-05 18:19:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t15 = qJD(1) + qJD(2);
	t16 = qJ(1) + qJ(2);
	t21 = sin(t16) * t15;
	t20 = cos(t16) * t15;
	t19 = r_i_i_C(1) * t21 + r_i_i_C(2) * t20;
	t18 = pkin(1) * qJD(1);
	t17 = -r_i_i_C(1) * t20 + r_i_i_C(2) * t21;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t18 + t19, t19, 0, 0, 0; -cos(qJ(1)) * t18 + t17, t17, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:19:03
	% EndTime: 2019-12-05 18:19:03
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->9), mult. (28->12), div. (0->0), fcn. (14->6), ass. (0->9)
	t22 = qJ(1) + qJ(2);
	t19 = pkin(8) + t22;
	t21 = qJD(1) + qJD(2);
	t26 = sin(t19) * t21;
	t25 = pkin(1) * qJD(1);
	t18 = cos(t19);
	t24 = r_i_i_C(1) * t26 + (sin(t22) * pkin(2) + r_i_i_C(2) * t18) * t21;
	t23 = r_i_i_C(2) * t26 + (-pkin(2) * cos(t22) - r_i_i_C(1) * t18) * t21;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t25 + t24, t24, 0, 0, 0; -cos(qJ(1)) * t25 + t23, t23, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:19:03
	% EndTime: 2019-12-05 18:19:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (109->16), mult. (66->19), div. (0->0), fcn. (40->8), ass. (0->14)
	t75 = r_i_i_C(1) * cos(pkin(9)) + pkin(3);
	t73 = r_i_i_C(2) * sin(pkin(9));
	t64 = qJ(1) + qJ(2);
	t61 = pkin(8) + t64;
	t59 = sin(t61);
	t63 = qJD(1) + qJD(2);
	t72 = t63 * t59;
	t60 = cos(t61);
	t71 = t63 * t60;
	t70 = -r_i_i_C(3) - qJ(4);
	t69 = pkin(1) * qJD(1);
	t68 = -t59 * qJD(4) + t75 * t72 + (pkin(2) * sin(t64) - t59 * t73 + t70 * t60) * t63;
	t67 = t71 * t73 + t60 * qJD(4) + (-pkin(2) * cos(t64) - t75 * t60 + t70 * t59) * t63;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t69 + t68, t68, 0, -t72, 0; -cos(qJ(1)) * t69 + t67, t67, 0, t71, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:19:04
	% EndTime: 2019-12-05 18:19:04
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (194->30), mult. (122->45), div. (0->0), fcn. (78->9), ass. (0->23)
	t105 = pkin(9) + qJ(5);
	t102 = cos(t105);
	t117 = qJD(5) * t102;
	t101 = sin(t105);
	t106 = qJD(1) + qJD(2);
	t121 = t101 * t106;
	t107 = qJ(1) + qJ(2);
	t103 = pkin(8) + t107;
	t98 = sin(t103);
	t99 = cos(t103);
	t125 = t98 * t117 + t99 * t121;
	t118 = qJD(5) * t101;
	t120 = t102 * t106;
	t124 = t99 * t118 + t98 * t120;
	t123 = pkin(1) * qJD(1);
	t122 = t106 * t98;
	t119 = t106 * (-pkin(7) - qJ(4));
	t114 = t98 * t118;
	t111 = t99 * t117;
	t100 = cos(pkin(9)) * pkin(4) + pkin(3);
	t110 = r_i_i_C(2) * t111 - t98 * qJD(4) + t100 * t122 + t99 * t119 + (-r_i_i_C(2) * t101 * t98 + pkin(2) * sin(t107) - r_i_i_C(3) * t99) * t106 + t124 * r_i_i_C(1);
	t109 = r_i_i_C(1) * t114 + t98 * t119 + t99 * qJD(4) + (-pkin(2) * cos(t107) - r_i_i_C(3) * t98 + (-r_i_i_C(1) * t102 - t100) * t99) * t106 + t125 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t101 - r_i_i_C(2) * t102) * qJD(5); sin(qJ(1)) * t123 + t110, t110, 0, -t122, (t99 * t120 - t114) * r_i_i_C(2) + t125 * r_i_i_C(1); -cos(qJ(1)) * t123 + t109, t109, 0, t106 * t99, t124 * r_i_i_C(2) + (t98 * t121 - t111) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end