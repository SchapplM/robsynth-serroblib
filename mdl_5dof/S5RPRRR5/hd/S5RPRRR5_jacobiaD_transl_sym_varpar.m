% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
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
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (sin(qJ(1)) * pkin(1) + r_i_i_C(1) * t7 + r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t8 + r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(9);
	t16 = qJ(3) + t18;
	t17 = qJD(1) + qJD(3);
	t22 = sin(t16) * t17;
	t21 = cos(t16) * t17;
	t20 = r_i_i_C(1) * t22 + r_i_i_C(2) * t21;
	t19 = -r_i_i_C(1) * t21 + r_i_i_C(2) * t22;
	t1 = [0, 0, 0, 0, 0; (sin(t18) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0; (-cos(t18) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t19, 0, t19, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (133->23), mult. (108->37), div. (0->0), fcn. (66->8), ass. (0->19)
	t92 = cos(qJ(4));
	t101 = qJD(4) * t92;
	t89 = qJD(1) + qJD(3);
	t91 = sin(qJ(4));
	t104 = t89 * t91;
	t90 = qJ(1) + pkin(9);
	t88 = qJ(3) + t90;
	t86 = sin(t88);
	t87 = cos(t88);
	t107 = t86 * t101 + t87 * t104;
	t102 = qJD(4) * t91;
	t103 = t89 * t92;
	t106 = t87 * t102 + t86 * t103;
	t105 = -pkin(7) - r_i_i_C(3);
	t98 = t86 * t102;
	t95 = t87 * t101;
	t94 = r_i_i_C(2) * t95 + (t105 * t87 + (-r_i_i_C(2) * t91 + pkin(3)) * t86) * t89 + t106 * r_i_i_C(1);
	t93 = r_i_i_C(1) * t98 + ((-r_i_i_C(1) * t92 - pkin(3)) * t87 + t105 * t86) * t89 + t107 * r_i_i_C(2);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t91 - r_i_i_C(2) * t92) * qJD(4), 0; (sin(t90) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t94, 0, t94, (t87 * t103 - t98) * r_i_i_C(2) + t107 * r_i_i_C(1), 0; (-cos(t90) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t93, 0, t93, t106 * r_i_i_C(2) + (t86 * t104 - t95) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:17:10
	% EndTime: 2019-12-05 18:17:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (251->32), mult. (166->47), div. (0->0), fcn. (105->10), ass. (0->30)
	t125 = qJ(1) + pkin(9);
	t120 = qJ(3) + t125;
	t117 = sin(t120);
	t118 = cos(t120);
	t126 = qJ(4) + qJ(5);
	t122 = cos(t126);
	t123 = qJD(4) + qJD(5);
	t146 = t122 * t123;
	t121 = sin(t126);
	t124 = qJD(1) + qJD(3);
	t147 = t121 * t124;
	t153 = -t117 * t147 + t118 * t146;
	t145 = t122 * t124;
	t148 = t121 * t123;
	t151 = t117 * t145 + t118 * t148;
	t150 = t117 * t146 + t118 * t147;
	t127 = sin(qJ(4));
	t141 = qJD(4) * t127 * pkin(4);
	t149 = t141 + (-r_i_i_C(3) - pkin(8) - pkin(7)) * t124;
	t144 = t124 * t127;
	t128 = cos(qJ(4));
	t142 = qJD(4) * t128;
	t140 = t117 * t148;
	t134 = (-r_i_i_C(1) * t121 - r_i_i_C(2) * t122) * t123;
	t133 = -t153 * r_i_i_C(1) + t151 * r_i_i_C(2);
	t132 = (t118 * t145 - t140) * r_i_i_C(2) + t150 * r_i_i_C(1);
	t119 = t128 * pkin(4) + pkin(3);
	t131 = t124 * t117 * t119 + t151 * r_i_i_C(1) + t153 * r_i_i_C(2) + t149 * t118;
	t130 = r_i_i_C(1) * t140 + (-r_i_i_C(1) * t122 - t119) * t118 * t124 + t149 * t117 + t150 * r_i_i_C(2);
	t1 = [0, 0, 0, t134 - t141, t134; (sin(t125) * pkin(2) + sin(qJ(1)) * pkin(1)) * qJD(1) + t131, 0, t131, (t117 * t142 + t118 * t144) * pkin(4) + t132, t132; (-cos(t125) * pkin(2) - cos(qJ(1)) * pkin(1)) * qJD(1) + t130, 0, t130, (t117 * t144 - t118 * t142) * pkin(4) + t133, t133;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end