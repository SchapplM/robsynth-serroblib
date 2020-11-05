% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:20
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t90 = cos(qJ(1));
	t89 = sin(qJ(1));
	t1 = [t90, -t89, 0, 0; t89, t90, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t91 = sin(pkin(6));
	t94 = sin(qJ(1));
	t102 = t94 * t91;
	t93 = sin(qJ(2));
	t101 = t94 * t93;
	t95 = cos(qJ(2));
	t100 = t94 * t95;
	t96 = cos(qJ(1));
	t99 = t96 * t91;
	t98 = t96 * t93;
	t97 = t96 * t95;
	t92 = cos(pkin(6));
	t1 = [-t92 * t101 + t97, -t92 * t100 - t98, t102, t96 * pkin(1) + pkin(8) * t102 + 0; t92 * t98 + t100, t92 * t97 - t101, -t99, t94 * pkin(1) - pkin(8) * t99 + 0; t91 * t93, t91 * t95, t92, t92 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t107 = sin(qJ(2));
	t108 = sin(qJ(1));
	t115 = t108 * t107;
	t109 = cos(qJ(2));
	t114 = t108 * t109;
	t110 = cos(qJ(1));
	t113 = t110 * t107;
	t112 = t110 * t109;
	t111 = pkin(2) * t107 - qJ(3) * t109;
	t106 = cos(pkin(6));
	t105 = sin(pkin(6));
	t104 = t109 * pkin(2) + t107 * qJ(3) + pkin(1);
	t103 = t105 * pkin(8) - t111 * t106;
	t1 = [t108 * t105, t106 * t115 - t112, t106 * t114 + t113, t103 * t108 + t104 * t110 + 0; -t110 * t105, -t106 * t113 - t114, -t106 * t112 + t115, -t103 * t110 + t104 * t108 + 0; t106, -t105 * t107, -t105 * t109, t106 * pkin(8) + t111 * t105 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t119 = sin(pkin(6));
	t124 = cos(qJ(4));
	t136 = t119 * t124;
	t121 = sin(qJ(4));
	t135 = t121 * t119;
	t122 = sin(qJ(2));
	t123 = sin(qJ(1));
	t134 = t123 * t122;
	t125 = cos(qJ(2));
	t133 = t123 * t125;
	t132 = t124 * t125;
	t126 = cos(qJ(1));
	t131 = t126 * t122;
	t130 = t126 * t125;
	t128 = pkin(2) + pkin(9);
	t129 = qJ(3) * t125 - t122 * t128;
	t127 = pkin(3) + pkin(8);
	t120 = cos(pkin(6));
	t118 = t122 * qJ(3) + t128 * t125 + pkin(1);
	t117 = -t120 * t130 + t134;
	t116 = t119 * t127 + t129 * t120;
	t1 = [t123 * t136 + (t120 * t133 + t131) * t121, (t120 * t132 - t135) * t123 + t124 * t131, -t120 * t134 + t130, t116 * t123 + t118 * t126 + 0; t117 * t121 - t126 * t136, t117 * t124 + t126 * t135, t120 * t131 + t133, -t116 * t126 + t118 * t123 + 0; t120 * t124 - t125 * t135, -t119 * t132 - t120 * t121, t119 * t122, -t129 * t119 + t127 * t120 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (65->36), mult. (112->59), div. (0->0), fcn. (146->10), ass. (0->29)
	t145 = sin(qJ(5));
	t147 = sin(qJ(2));
	t164 = t145 * t147;
	t143 = sin(pkin(6));
	t146 = sin(qJ(4));
	t163 = t146 * t143;
	t151 = cos(qJ(2));
	t162 = t146 * t151;
	t149 = cos(qJ(5));
	t161 = t147 * t149;
	t152 = cos(qJ(1));
	t160 = t147 * t152;
	t148 = sin(qJ(1));
	t159 = t148 * t147;
	t150 = cos(qJ(4));
	t158 = t150 * t151;
	t157 = t152 * t151;
	t142 = t146 * pkin(4) - pkin(10) * t150 + qJ(3);
	t153 = pkin(2) + pkin(9);
	t156 = t142 * t151 - t153 * t147;
	t155 = t150 * pkin(4) + t146 * pkin(10) + pkin(3) + pkin(8);
	t144 = cos(pkin(6));
	t139 = t150 * t143 + t144 * t162;
	t154 = t139 * t145 + t144 * t161;
	t141 = -t146 * t164 + t149 * t151;
	t140 = -t143 * t162 + t144 * t150;
	t138 = t142 * t147 + t153 * t151 + pkin(1);
	t137 = t143 * t155 + t156 * t144;
	t1 = [(t139 * t148 + t146 * t160) * t149 + (-t144 * t159 + t157) * t145, t152 * t141 - t154 * t148, (-t144 * t158 + t163) * t148 - t150 * t160, t137 * t148 + t138 * t152 + 0; (-t139 * t149 + t144 * t164) * t152 + t148 * (t145 * t151 + t146 * t161), t148 * t141 + t154 * t152, (t144 * t157 - t159) * t150 - t152 * t163, -t137 * t152 + t138 * t148 + 0; t140 * t149 + t143 * t164, -t140 * t145 + t143 * t161, t143 * t158 + t144 * t146, -t156 * t143 + t155 * t144 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:20:47
	% EndTime: 2020-11-04 22:20:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (98->38), mult. (125->56), div. (0->0), fcn. (159->12), ass. (0->33)
	t177 = sin(pkin(6));
	t180 = sin(qJ(2));
	t197 = t177 * t180;
	t179 = sin(qJ(4));
	t196 = t179 * t177;
	t183 = cos(qJ(2));
	t195 = t179 * t183;
	t184 = cos(qJ(1));
	t194 = t180 * t184;
	t181 = sin(qJ(1));
	t193 = t181 * t180;
	t182 = cos(qJ(4));
	t192 = t182 * t183;
	t191 = t184 * t183;
	t172 = sin(qJ(5)) * pkin(5) + pkin(2) + pkin(9);
	t173 = cos(qJ(5)) * pkin(5) + pkin(4);
	t185 = pkin(11) + pkin(10);
	t187 = t173 * t179 - t185 * t182 + qJ(3);
	t190 = t172 * t180 - t183 * t187;
	t178 = cos(pkin(6));
	t168 = t182 * t177 + t178 * t195;
	t189 = t168 * t181 + t179 * t194;
	t188 = -t168 * t184 + t179 * t193;
	t186 = t173 * t182 + t185 * t179 + pkin(3) + pkin(8);
	t176 = qJ(5) + qJ(6);
	t175 = cos(t176);
	t174 = sin(t176);
	t171 = -t178 * t193 + t191;
	t170 = t178 * t194 + t181 * t183;
	t169 = -t177 * t195 + t178 * t182;
	t166 = t172 * t183 + t187 * t180 + pkin(1);
	t165 = -t177 * t186 + t190 * t178;
	t1 = [t171 * t174 + t189 * t175, t171 * t175 - t189 * t174, (-t178 * t192 + t196) * t181 - t182 * t194, -t165 * t181 + t166 * t184 + 0; t170 * t174 + t188 * t175, t170 * t175 - t188 * t174, (t178 * t191 - t193) * t182 - t184 * t196, t165 * t184 + t166 * t181 + 0; t169 * t175 + t174 * t197, -t169 * t174 + t175 * t197, t177 * t192 + t178 * t179, t190 * t177 + t186 * t178 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end