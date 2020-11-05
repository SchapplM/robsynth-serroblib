% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t88 = cos(qJ(1));
	t87 = sin(qJ(1));
	t1 = [t88, -t87, 0, 0; t87, t88, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t89 = sin(pkin(6));
	t92 = sin(qJ(1));
	t100 = t92 * t89;
	t91 = sin(qJ(2));
	t99 = t92 * t91;
	t93 = cos(qJ(2));
	t98 = t92 * t93;
	t94 = cos(qJ(1));
	t97 = t94 * t89;
	t96 = t94 * t91;
	t95 = t94 * t93;
	t90 = cos(pkin(6));
	t1 = [-t90 * t99 + t95, -t90 * t98 - t96, t100, pkin(1) * t94 + pkin(8) * t100 + 0; t90 * t96 + t98, t90 * t95 - t99, -t97, pkin(1) * t92 - pkin(8) * t97 + 0; t89 * t91, t89 * t93, t90, pkin(8) * t90 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t105 = sin(qJ(2));
	t106 = sin(qJ(1));
	t113 = t106 * t105;
	t107 = cos(qJ(2));
	t112 = t106 * t107;
	t108 = cos(qJ(1));
	t111 = t108 * t105;
	t110 = t108 * t107;
	t109 = pkin(2) * t105 - qJ(3) * t107;
	t104 = cos(pkin(6));
	t103 = sin(pkin(6));
	t102 = t107 * pkin(2) + t105 * qJ(3) + pkin(1);
	t101 = t103 * pkin(8) - t109 * t104;
	t1 = [t106 * t103, t104 * t113 - t110, t104 * t112 + t111, t101 * t106 + t102 * t108 + 0; -t108 * t103, -t104 * t111 - t112, -t104 * t110 + t113, -t101 * t108 + t102 * t106 + 0; t104, -t103 * t105, -t103 * t107, t104 * pkin(8) + t109 * t103 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t117 = sin(pkin(6));
	t122 = cos(qJ(4));
	t134 = t117 * t122;
	t119 = sin(qJ(4));
	t133 = t119 * t117;
	t120 = sin(qJ(2));
	t121 = sin(qJ(1));
	t132 = t121 * t120;
	t123 = cos(qJ(2));
	t131 = t121 * t123;
	t130 = t122 * t123;
	t124 = cos(qJ(1));
	t129 = t124 * t120;
	t128 = t124 * t123;
	t126 = pkin(2) + pkin(9);
	t127 = qJ(3) * t123 - t120 * t126;
	t125 = pkin(3) + pkin(8);
	t118 = cos(pkin(6));
	t116 = t120 * qJ(3) + t126 * t123 + pkin(1);
	t115 = -t118 * t128 + t132;
	t114 = t117 * t125 + t118 * t127;
	t1 = [t121 * t134 + (t118 * t131 + t129) * t119, (t118 * t130 - t133) * t121 + t122 * t129, -t118 * t132 + t128, t114 * t121 + t116 * t124 + 0; t115 * t119 - t124 * t134, t115 * t122 + t124 * t133, t118 * t129 + t131, -t114 * t124 + t116 * t121 + 0; t118 * t122 - t123 * t133, -t117 * t130 - t118 * t119, t117 * t120, -t127 * t117 + t125 * t118 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (65->36), mult. (112->59), div. (0->0), fcn. (146->10), ass. (0->29)
	t143 = sin(qJ(5));
	t145 = sin(qJ(2));
	t162 = t143 * t145;
	t141 = sin(pkin(6));
	t144 = sin(qJ(4));
	t161 = t144 * t141;
	t149 = cos(qJ(2));
	t160 = t144 * t149;
	t147 = cos(qJ(5));
	t159 = t145 * t147;
	t150 = cos(qJ(1));
	t158 = t145 * t150;
	t146 = sin(qJ(1));
	t157 = t146 * t145;
	t148 = cos(qJ(4));
	t156 = t148 * t149;
	t155 = t150 * t149;
	t140 = t144 * pkin(4) - pkin(10) * t148 + qJ(3);
	t151 = pkin(2) + pkin(9);
	t154 = t140 * t149 - t151 * t145;
	t153 = t148 * pkin(4) + t144 * pkin(10) + pkin(3) + pkin(8);
	t142 = cos(pkin(6));
	t137 = t148 * t141 + t142 * t160;
	t152 = t137 * t143 + t142 * t159;
	t139 = -t144 * t162 + t149 * t147;
	t138 = -t141 * t160 + t142 * t148;
	t136 = t140 * t145 + t151 * t149 + pkin(1);
	t135 = t141 * t153 + t154 * t142;
	t1 = [(t137 * t146 + t144 * t158) * t147 + (-t142 * t157 + t155) * t143, t150 * t139 - t152 * t146, (-t142 * t156 + t161) * t146 - t148 * t158, t135 * t146 + t136 * t150 + 0; (-t137 * t147 + t142 * t162) * t150 + t146 * (t149 * t143 + t144 * t159), t146 * t139 + t152 * t150, (t142 * t155 - t157) * t148 - t150 * t161, -t135 * t150 + t136 * t146 + 0; t138 * t147 + t141 * t162, -t138 * t143 + t141 * t159, t141 * t156 + t142 * t144, -t154 * t141 + t153 * t142 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:16:48
	% EndTime: 2020-11-04 22:16:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (86->39), mult. (125->61), div. (0->0), fcn. (159->10), ass. (0->31)
	t174 = sin(qJ(5));
	t176 = sin(qJ(2));
	t192 = t174 * t176;
	t171 = sin(pkin(6));
	t175 = sin(qJ(4));
	t191 = t175 * t171;
	t180 = cos(qJ(2));
	t190 = t175 * t180;
	t178 = cos(qJ(5));
	t189 = t176 * t178;
	t181 = cos(qJ(1));
	t188 = t176 * t181;
	t177 = sin(qJ(1));
	t187 = t177 * t176;
	t179 = cos(qJ(4));
	t186 = t179 * t180;
	t185 = t181 * t180;
	t170 = t178 * pkin(5) + pkin(4);
	t173 = qJ(6) + pkin(10);
	t165 = t170 * t175 - t173 * t179 + qJ(3);
	t169 = t174 * pkin(5) + pkin(2) + pkin(9);
	t184 = t165 * t180 - t169 * t176;
	t172 = cos(pkin(6));
	t166 = t179 * t171 + t172 * t190;
	t183 = t166 * t174 + t172 * t189;
	t182 = t170 * t179 + t173 * t175 + pkin(3) + pkin(8);
	t168 = -t175 * t192 + t180 * t178;
	t167 = -t171 * t190 + t172 * t179;
	t164 = t165 * t176 + t169 * t180 + pkin(1);
	t163 = t171 * t182 + t184 * t172;
	t1 = [(t166 * t177 + t175 * t188) * t178 + (-t172 * t187 + t185) * t174, t181 * t168 - t183 * t177, (-t172 * t186 + t191) * t177 - t179 * t188, t163 * t177 + t164 * t181 + 0; (-t166 * t178 + t172 * t192) * t181 + t177 * (t180 * t174 + t175 * t189), t177 * t168 + t183 * t181, (t172 * t185 - t187) * t179 - t181 * t191, -t163 * t181 + t164 * t177 + 0; t167 * t178 + t171 * t192, -t167 * t174 + t171 * t189, t171 * t186 + t172 * t175, -t184 * t171 + t182 * t172 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end