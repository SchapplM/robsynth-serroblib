% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
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
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
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
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t24 = qJD(1) + qJD(2);
	t33 = pkin(2) * t24;
	t21 = qJD(3) + t24;
	t25 = qJ(1) + qJ(2);
	t23 = qJ(3) + t25;
	t32 = sin(t23) * t21;
	t31 = cos(t23) * t21;
	t30 = r_i_i_C(1) * t32 + r_i_i_C(2) * t31;
	t29 = pkin(1) * qJD(1);
	t28 = sin(t25) * t33 + t30;
	t27 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t32;
	t26 = -cos(t25) * t33 + t27;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t29 + t28, t28, t30, 0, 0; -cos(qJ(1)) * t29 + t26, t26, t27, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (225->24), mult. (148->37), div. (0->0), fcn. (90->8), ass. (0->24)
	t98 = sin(qJ(4));
	t111 = qJD(4) * t98;
	t96 = qJD(1) + qJD(2);
	t93 = qJD(3) + t96;
	t99 = cos(qJ(4));
	t113 = t93 * t99;
	t97 = qJ(1) + qJ(2);
	t95 = qJ(3) + t97;
	t91 = sin(t95);
	t92 = cos(t95);
	t118 = t92 * t111 + t91 * t113;
	t110 = qJD(4) * t99;
	t114 = t93 * t98;
	t117 = t91 * t110 + t92 * t114;
	t116 = -pkin(8) - r_i_i_C(3);
	t115 = pkin(2) * t96;
	t112 = pkin(1) * qJD(1);
	t107 = t91 * t111;
	t104 = t92 * t110;
	t103 = r_i_i_C(2) * t104 + (t116 * t92 + (-r_i_i_C(2) * t98 + pkin(3)) * t91) * t93 + t118 * r_i_i_C(1);
	t102 = sin(t97) * t115 + t103;
	t101 = r_i_i_C(1) * t107 + ((-r_i_i_C(1) * t99 - pkin(3)) * t92 + t116 * t91) * t93 + t117 * r_i_i_C(2);
	t100 = -cos(t97) * t115 + t101;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t98 - r_i_i_C(2) * t99) * qJD(4), 0; sin(qJ(1)) * t112 + t102, t102, t103, (t92 * t113 - t107) * r_i_i_C(2) + t117 * r_i_i_C(1), 0; -cos(qJ(1)) * t112 + t100, t100, t101, t118 * r_i_i_C(2) + (t91 * t114 - t104) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (371->33), mult. (214->47), div. (0->0), fcn. (135->10), ass. (0->35)
	t133 = qJ(1) + qJ(2);
	t129 = qJ(3) + t133;
	t122 = sin(t129);
	t123 = cos(t129);
	t132 = qJ(4) + qJ(5);
	t127 = cos(t132);
	t130 = qJD(4) + qJD(5);
	t152 = t127 * t130;
	t131 = qJD(1) + qJD(2);
	t125 = qJD(3) + t131;
	t126 = sin(t132);
	t157 = t125 * t126;
	t164 = -t122 * t157 + t123 * t152;
	t162 = t122 * t152 + t123 * t157;
	t153 = t126 * t130;
	t156 = t125 * t127;
	t161 = t122 * t156 + t123 * t153;
	t134 = sin(qJ(4));
	t150 = qJD(4) * t134 * pkin(4);
	t160 = t150 + (-r_i_i_C(3) - pkin(9) - pkin(8)) * t125;
	t159 = pkin(2) * t131;
	t158 = pkin(1) * qJD(1);
	t155 = t125 * t134;
	t135 = cos(qJ(4));
	t151 = qJD(4) * t135;
	t148 = t122 * t153;
	t143 = (-r_i_i_C(1) * t126 - r_i_i_C(2) * t127) * t130;
	t142 = -t164 * r_i_i_C(1) + t161 * r_i_i_C(2);
	t141 = (t123 * t156 - t148) * r_i_i_C(2) + t162 * r_i_i_C(1);
	t124 = t135 * pkin(4) + pkin(3);
	t140 = t125 * t122 * t124 + t161 * r_i_i_C(1) + t164 * r_i_i_C(2) + t160 * t123;
	t139 = sin(t133) * t159 + t140;
	t138 = r_i_i_C(1) * t148 + (-r_i_i_C(1) * t127 - t124) * t123 * t125 + t160 * t122 + t162 * r_i_i_C(2);
	t137 = -cos(t133) * t159 + t138;
	t1 = [0, 0, 0, t143 - t150, t143; sin(qJ(1)) * t158 + t139, t139, t140, (t122 * t151 + t123 * t155) * pkin(4) + t141, t141; -cos(qJ(1)) * t158 + t137, t137, t138, (t122 * t155 - t123 * t151) * pkin(4) + t142, t142;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end