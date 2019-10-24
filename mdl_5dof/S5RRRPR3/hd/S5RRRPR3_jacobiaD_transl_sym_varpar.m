% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:50
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
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
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (159->23), mult. (140->28), div. (0->0), fcn. (91->8), ass. (0->22)
	t105 = qJ(3) + pkin(9);
	t100 = sin(t105);
	t101 = cos(t105);
	t114 = r_i_i_C(1) * t100 + r_i_i_C(2) * t101 + sin(qJ(3)) * pkin(3);
	t134 = t114 * qJD(3);
	t133 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t101;
	t132 = pkin(2) + t133;
	t104 = qJD(1) + qJD(2);
	t125 = r_i_i_C(2) * t100;
	t131 = -t104 * t125 - qJD(4);
	t130 = t104 * (-qJ(4) - pkin(7)) + t134;
	t129 = qJD(3) * (-t125 + t133);
	t122 = pkin(1) * qJD(1);
	t106 = qJ(1) + qJ(2);
	t102 = sin(t106);
	t121 = t104 * t102;
	t103 = cos(t106);
	t120 = t104 * t103;
	t112 = t114 * t104;
	t111 = (-t132 * t104 - t131) * t103 + (-r_i_i_C(3) * t104 + t130) * t102;
	t110 = -r_i_i_C(3) * t120 + t131 * t102 + t130 * t103 + t132 * t121;
	t1 = [0, 0, -t134, 0, 0; sin(qJ(1)) * t122 + t110, t110, t102 * t129 + t103 * t112, -t121, 0; -cos(qJ(1)) * t122 + t111, t111, t102 * t112 - t103 * t129, t120, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:56
	% EndTime: 2019-10-24 10:50:56
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (274->40), mult. (188->52), div. (0->0), fcn. (122->10), ass. (0->34)
	t131 = qJ(3) + pkin(9);
	t125 = qJ(5) + t131;
	t121 = sin(t125);
	t122 = cos(t125);
	t132 = qJ(1) + qJ(2);
	t127 = cos(t132);
	t130 = qJD(1) + qJD(2);
	t147 = t130 * t127;
	t126 = sin(t132);
	t129 = qJD(3) + qJD(5);
	t150 = t126 * t129;
	t157 = t121 * t147 + t122 * t150;
	t148 = t130 * t126;
	t149 = t127 * t129;
	t156 = t121 * t149 + t122 * t148;
	t155 = r_i_i_C(1) * t121;
	t154 = r_i_i_C(1) * t122;
	t153 = r_i_i_C(2) * t121;
	t152 = r_i_i_C(2) * t122;
	t151 = pkin(1) * qJD(1);
	t146 = t121 * t150;
	t141 = t122 * t149;
	t140 = t157 * r_i_i_C(1) + t147 * t152;
	t139 = t156 * r_i_i_C(2) + t148 * t155;
	t119 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t131);
	t138 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t131);
	t137 = (-t152 - t155) * t129;
	t114 = t119 * qJD(3);
	t118 = pkin(2) - t138;
	t128 = -pkin(8) - qJ(4) - pkin(7);
	t136 = t118 * t148 + r_i_i_C(2) * t141 + t128 * t147 + (-r_i_i_C(3) * t130 - t114) * t127 + (-t130 * t153 - qJD(4)) * t126 + t156 * r_i_i_C(1);
	t135 = -t126 * t114 + r_i_i_C(1) * t146 + t128 * t148 + t127 * qJD(4) + (-r_i_i_C(3) * t126 + (-t118 - t154) * t127) * t130 + t157 * r_i_i_C(2);
	t115 = t138 * qJD(3);
	t1 = [0, 0, t137 + t114, 0, t137; sin(qJ(1)) * t151 + t136, t136, -t119 * t147 + (-t129 * t153 - t115) * t126 + t140, -t148, -r_i_i_C(2) * t146 + t140; -cos(qJ(1)) * t151 + t135, t135, -t119 * t148 + (-t129 * t154 + t115) * t127 + t139, t147, -r_i_i_C(1) * t141 + t139;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end