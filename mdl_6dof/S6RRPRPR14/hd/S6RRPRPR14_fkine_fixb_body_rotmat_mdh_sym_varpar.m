% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t87 = cos(qJ(1));
	t86 = sin(qJ(1));
	t1 = [t87, -t86, 0, 0; t86, t87, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t88 = sin(pkin(6));
	t91 = sin(qJ(1));
	t99 = t91 * t88;
	t90 = sin(qJ(2));
	t98 = t91 * t90;
	t92 = cos(qJ(2));
	t97 = t91 * t92;
	t93 = cos(qJ(1));
	t96 = t93 * t88;
	t95 = t93 * t90;
	t94 = t93 * t92;
	t89 = cos(pkin(6));
	t1 = [-t89 * t98 + t94, -t89 * t97 - t95, t99, pkin(1) * t93 + pkin(8) * t99 + 0; t89 * t95 + t97, t89 * t94 - t98, -t96, pkin(1) * t91 - pkin(8) * t96 + 0; t88 * t90, t88 * t92, t89, pkin(8) * t89 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t104 = sin(qJ(2));
	t105 = sin(qJ(1));
	t112 = t105 * t104;
	t106 = cos(qJ(2));
	t111 = t105 * t106;
	t107 = cos(qJ(1));
	t110 = t107 * t104;
	t109 = t107 * t106;
	t108 = pkin(2) * t104 - qJ(3) * t106;
	t103 = cos(pkin(6));
	t102 = sin(pkin(6));
	t101 = t106 * pkin(2) + t104 * qJ(3) + pkin(1);
	t100 = t102 * pkin(8) - t108 * t103;
	t1 = [t105 * t102, t103 * t112 - t109, t103 * t111 + t110, t100 * t105 + t101 * t107 + 0; -t107 * t102, -t103 * t110 - t111, -t103 * t109 + t112, -t100 * t107 + t101 * t105 + 0; t103, -t102 * t104, -t102 * t106, t103 * pkin(8) + t108 * t102 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t116 = sin(pkin(6));
	t121 = cos(qJ(4));
	t133 = t116 * t121;
	t118 = sin(qJ(4));
	t132 = t118 * t116;
	t119 = sin(qJ(2));
	t120 = sin(qJ(1));
	t131 = t120 * t119;
	t122 = cos(qJ(2));
	t130 = t120 * t122;
	t129 = t121 * t122;
	t123 = cos(qJ(1));
	t128 = t123 * t119;
	t127 = t123 * t122;
	t125 = pkin(2) + pkin(9);
	t126 = qJ(3) * t122 - t119 * t125;
	t124 = pkin(3) + pkin(8);
	t117 = cos(pkin(6));
	t115 = t119 * qJ(3) + t125 * t122 + pkin(1);
	t114 = -t117 * t127 + t131;
	t113 = t116 * t124 + t117 * t126;
	t1 = [t120 * t133 + (t117 * t130 + t128) * t118, (t117 * t129 - t132) * t120 + t121 * t128, -t117 * t131 + t127, t113 * t120 + t115 * t123 + 0; t114 * t118 - t123 * t133, t114 * t121 + t123 * t132, t117 * t128 + t130, -t113 * t123 + t115 * t120 + 0; t117 * t121 - t122 * t132, -t116 * t129 - t117 * t118, t116 * t119, -t116 * t126 + t124 * t117 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (52->28), mult. (78->41), div. (0->0), fcn. (99->8), ass. (0->23)
	t138 = sin(pkin(6));
	t143 = cos(qJ(4));
	t155 = t138 * t143;
	t140 = sin(qJ(4));
	t154 = t140 * t138;
	t141 = sin(qJ(2));
	t142 = sin(qJ(1));
	t153 = t142 * t141;
	t144 = cos(qJ(2));
	t152 = t142 * t144;
	t151 = t143 * t144;
	t145 = cos(qJ(1));
	t150 = t145 * t141;
	t149 = t145 * t144;
	t137 = t140 * pkin(4) - qJ(5) * t143 + qJ(3);
	t146 = pkin(2) + pkin(9);
	t148 = t137 * t144 - t146 * t141;
	t147 = pkin(4) * t143 + qJ(5) * t140 + pkin(3) + pkin(8);
	t139 = cos(pkin(6));
	t136 = t139 * t149 - t153;
	t135 = t137 * t141 + t146 * t144 + pkin(1);
	t134 = t138 * t147 + t148 * t139;
	t1 = [-t139 * t153 + t149, (-t139 * t152 - t150) * t140 - t142 * t155, (-t139 * t151 + t154) * t142 - t143 * t150, t134 * t142 + t135 * t145 + 0; t139 * t150 + t152, t136 * t140 + t145 * t155, t136 * t143 - t145 * t154, -t134 * t145 + t135 * t142 + 0; t138 * t141, -t139 * t143 + t144 * t154, t138 * t151 + t139 * t140, -t148 * t138 + t147 * t139 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:46
	% EndTime: 2020-11-04 22:12:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (78->37), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->31)
	t164 = sin(pkin(6));
	t171 = cos(qJ(4));
	t186 = t164 * t171;
	t166 = sin(qJ(6));
	t168 = sin(qJ(2));
	t185 = t166 * t168;
	t167 = sin(qJ(4));
	t184 = t167 * t164;
	t170 = cos(qJ(6));
	t183 = t168 * t170;
	t182 = t168 * t171;
	t169 = sin(qJ(1));
	t181 = t169 * t168;
	t172 = cos(qJ(2));
	t180 = t171 * t172;
	t173 = cos(qJ(1));
	t179 = t173 * t168;
	t178 = t173 * t172;
	t174 = pkin(4) + pkin(10);
	t161 = -qJ(5) * t171 + t174 * t167 + qJ(3);
	t163 = pkin(2) + pkin(5) + pkin(9);
	t177 = t161 * t172 - t163 * t168;
	t165 = cos(pkin(6));
	t158 = t165 * t180 - t184;
	t176 = t158 * t166 + t165 * t183;
	t175 = qJ(5) * t167 + t174 * t171 + pkin(3) + pkin(8);
	t160 = -t166 * t182 + t170 * t172;
	t159 = t164 * t180 + t165 * t167;
	t157 = t161 * t168 + t163 * t172 + pkin(1);
	t156 = t164 * t175 + t177 * t165;
	t1 = [t173 * t160 - t176 * t169, -(-t165 * t181 + t178) * t166 + (-t158 * t169 - t171 * t179) * t170, t169 * t186 + (t169 * t165 * t172 + t179) * t167, t156 * t169 + t157 * t173 + 0; t169 * t160 + t176 * t173, (t158 * t170 - t165 * t185) * t173 - t169 * (t166 * t172 + t170 * t182), -t173 * t186 + (-t165 * t178 + t181) * t167, -t156 * t173 + t157 * t169 + 0; t159 * t166 + t164 * t183, t159 * t170 - t164 * t185, t165 * t171 - t172 * t184, -t177 * t164 + t175 * t165 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end