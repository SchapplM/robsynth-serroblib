% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:51
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:40
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
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:40
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
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (93->20), mult. (104->34), div. (0->0), fcn. (64->6), ass. (0->19)
	t89 = sin(qJ(3));
	t100 = qJD(3) * t89;
	t87 = qJD(1) + qJD(2);
	t90 = cos(qJ(3));
	t102 = t87 * t90;
	t88 = qJ(1) + qJ(2);
	t85 = sin(t88);
	t86 = cos(t88);
	t106 = t86 * t100 + t85 * t102;
	t103 = t87 * t89;
	t99 = qJD(3) * t90;
	t105 = t86 * t103 + t85 * t99;
	t104 = -pkin(7) - r_i_i_C(3);
	t101 = pkin(1) * qJD(1);
	t96 = t85 * t100;
	t93 = t86 * t99;
	t92 = r_i_i_C(2) * t93 + (t104 * t86 + (-r_i_i_C(2) * t89 + pkin(2)) * t85) * t87 + t106 * r_i_i_C(1);
	t91 = r_i_i_C(1) * t96 + ((-r_i_i_C(1) * t90 - pkin(2)) * t86 + t104 * t85) * t87 + t105 * r_i_i_C(2);
	t1 = [0, 0, (-r_i_i_C(1) * t89 - r_i_i_C(2) * t90) * qJD(3), 0, 0; sin(qJ(1)) * t101 + t92, t92, (t86 * t102 - t96) * r_i_i_C(2) + t105 * r_i_i_C(1), 0, 0; -cos(qJ(1)) * t101 + t91, t91, t106 * r_i_i_C(2) + (t85 * t103 - t93) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (195->30), mult. (162->43), div. (0->0), fcn. (103->8), ass. (0->32)
	t124 = qJ(1) + qJ(2);
	t118 = sin(t124);
	t122 = qJD(1) + qJD(2);
	t145 = t118 * t122;
	t120 = cos(t124);
	t143 = t120 * t122;
	t123 = qJ(3) + qJ(4);
	t117 = sin(t123);
	t119 = cos(t123);
	t121 = qJD(3) + qJD(4);
	t146 = t118 * t121;
	t152 = t117 * t143 + t119 * t146;
	t144 = t120 * t121;
	t151 = t117 * t144 + t119 * t145;
	t125 = sin(qJ(3));
	t139 = qJD(3) * t125 * pkin(3);
	t150 = t139 + (-r_i_i_C(3) - pkin(8) - pkin(7)) * t122;
	t149 = r_i_i_C(1) * t117;
	t148 = r_i_i_C(2) * t119;
	t147 = pkin(1) * qJD(1);
	t142 = t122 * t125;
	t126 = cos(qJ(3));
	t140 = qJD(3) * t126;
	t138 = t117 * t146;
	t133 = t119 * t144;
	t132 = (-t148 - t149) * t121;
	t131 = -r_i_i_C(1) * t133 + r_i_i_C(2) * t151 + t145 * t149;
	t130 = r_i_i_C(1) * t152 - r_i_i_C(2) * t138 + t143 * t148;
	t116 = t126 * pkin(3) + pkin(2);
	t129 = t116 * t145 + (-t117 * t145 + t133) * r_i_i_C(2) + t150 * t120 + t151 * r_i_i_C(1);
	t128 = r_i_i_C(1) * t138 + (-r_i_i_C(1) * t119 - t116) * t143 + t150 * t118 + t152 * r_i_i_C(2);
	t1 = [0, 0, t132 - t139, t132, 0; sin(qJ(1)) * t147 + t129, t129, (t118 * t140 + t120 * t142) * pkin(3) + t130, t130, 0; -cos(qJ(1)) * t147 + t128, t128, (t118 * t142 - t120 * t140) * pkin(3) + t131, t131, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:51:40
	% EndTime: 2019-10-24 10:51:41
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (264->45), mult. (202->59), div. (0->0), fcn. (131->8), ass. (0->35)
	t127 = qJ(3) + qJ(4);
	t120 = sin(t127);
	t128 = qJ(1) + qJ(2);
	t121 = sin(t128);
	t123 = cos(t128);
	t126 = qJD(1) + qJD(2);
	t144 = t126 * t123;
	t122 = cos(t127);
	t125 = qJD(3) + qJD(4);
	t147 = t122 * t125;
	t153 = t120 * t144 + t121 * t147;
	t145 = t126 * t121;
	t146 = t123 * t125;
	t152 = t120 * t146 + t122 * t145;
	t151 = r_i_i_C(2) * t122;
	t150 = pkin(1) * qJD(1);
	t149 = pkin(3) * qJD(3);
	t148 = t120 * t125;
	t129 = sin(qJ(3));
	t143 = t129 * t149;
	t142 = t121 * t148;
	t141 = t120 * t145;
	t136 = t122 * t146;
	t135 = t153 * r_i_i_C(1) + t144 * t151;
	t134 = r_i_i_C(1) * t141 + t152 * r_i_i_C(2);
	t133 = (-t151 + (-pkin(4) - r_i_i_C(1)) * t120) * t125;
	t101 = -pkin(4) * t148 - t143;
	t130 = cos(qJ(3));
	t117 = t130 * pkin(3) + pkin(4) * t122 + pkin(2);
	t124 = -qJ(5) - pkin(8) - pkin(7);
	t132 = t117 * t145 + r_i_i_C(2) * t136 + t124 * t144 + (-r_i_i_C(3) * t126 - t101) * t123 + (-r_i_i_C(2) * t120 * t126 - qJD(5)) * t121 + t152 * r_i_i_C(1);
	t131 = -t121 * t101 + r_i_i_C(1) * t142 + t124 * t145 + t123 * qJD(5) + (-r_i_i_C(3) * t121 + (-r_i_i_C(1) * t122 - t117) * t123) * t126 + t153 * r_i_i_C(2);
	t118 = -t129 * pkin(3) - pkin(4) * t120;
	t102 = -pkin(4) * t147 - t130 * t149;
	t1 = [0, 0, t133 - t143, t133, 0; sin(qJ(1)) * t150 + t132, t132, -t118 * t144 + (-r_i_i_C(2) * t148 - t102) * t121 + t135, t153 * pkin(4) - r_i_i_C(2) * t142 + t135, -t145; -cos(qJ(1)) * t150 + t131, t131, -t118 * t145 + (-r_i_i_C(1) * t147 + t102) * t123 + t134, -r_i_i_C(1) * t136 + (-t136 + t141) * pkin(4) + t134, t144;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end