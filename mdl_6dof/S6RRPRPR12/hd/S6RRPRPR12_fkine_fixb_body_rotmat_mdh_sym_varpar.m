% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t1 = [t83, -t82, 0, 0; t82, t83, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t84 = sin(pkin(6));
	t87 = sin(qJ(1));
	t95 = t87 * t84;
	t86 = sin(qJ(2));
	t94 = t87 * t86;
	t88 = cos(qJ(2));
	t93 = t87 * t88;
	t89 = cos(qJ(1));
	t92 = t89 * t84;
	t91 = t89 * t86;
	t90 = t89 * t88;
	t85 = cos(pkin(6));
	t1 = [-t85 * t94 + t90, -t85 * t93 - t91, t95, t89 * pkin(1) + pkin(8) * t95 + 0; t85 * t91 + t93, t85 * t90 - t94, -t92, t87 * pkin(1) - pkin(8) * t92 + 0; t84 * t86, t84 * t88, t85, t85 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t100 = sin(qJ(2));
	t101 = sin(qJ(1));
	t108 = t101 * t100;
	t102 = cos(qJ(2));
	t107 = t101 * t102;
	t103 = cos(qJ(1));
	t106 = t103 * t100;
	t105 = t103 * t102;
	t104 = pkin(2) * t100 - qJ(3) * t102;
	t99 = cos(pkin(6));
	t98 = sin(pkin(6));
	t97 = t102 * pkin(2) + t100 * qJ(3) + pkin(1);
	t96 = t98 * pkin(8) - t104 * t99;
	t1 = [t101 * t98, t99 * t108 - t105, t99 * t107 + t106, t96 * t101 + t97 * t103 + 0; -t103 * t98, -t99 * t106 - t107, -t99 * t105 + t108, t97 * t101 - t96 * t103 + 0; t99, -t98 * t100, -t98 * t102, t99 * pkin(8) + t104 * t98 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t112 = sin(pkin(6));
	t117 = cos(qJ(4));
	t129 = t112 * t117;
	t114 = sin(qJ(4));
	t128 = t114 * t112;
	t115 = sin(qJ(2));
	t116 = sin(qJ(1));
	t127 = t116 * t115;
	t118 = cos(qJ(2));
	t126 = t116 * t118;
	t125 = t117 * t118;
	t119 = cos(qJ(1));
	t124 = t119 * t115;
	t123 = t119 * t118;
	t121 = pkin(2) + pkin(9);
	t122 = qJ(3) * t118 - t115 * t121;
	t120 = pkin(3) + pkin(8);
	t113 = cos(pkin(6));
	t111 = t115 * qJ(3) + t121 * t118 + pkin(1);
	t110 = -t113 * t123 + t127;
	t109 = t112 * t120 + t122 * t113;
	t1 = [t116 * t129 + (t113 * t126 + t124) * t114, (t113 * t125 - t128) * t116 + t117 * t124, -t113 * t127 + t123, t109 * t116 + t111 * t119 + 0; t110 * t114 - t119 * t129, t110 * t117 + t119 * t128, t113 * t124 + t126, -t109 * t119 + t111 * t116 + 0; t113 * t117 - t118 * t128, -t112 * t125 - t113 * t114, t112 * t115, -t122 * t112 + t120 * t113 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (61->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t139 = sin(pkin(6));
	t142 = sin(qJ(1));
	t154 = t139 * t142;
	t144 = cos(qJ(2));
	t153 = t139 * t144;
	t145 = cos(qJ(1));
	t152 = t139 * t145;
	t141 = sin(qJ(2));
	t151 = t142 * t141;
	t150 = t142 * t144;
	t149 = t145 * t141;
	t148 = t145 * t144;
	t147 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(8);
	t134 = sin(qJ(4)) * pkin(4) + qJ(3);
	t137 = qJ(5) + pkin(2) + pkin(9);
	t146 = t134 * t144 - t137 * t141;
	t140 = cos(pkin(6));
	t138 = qJ(4) + pkin(11);
	t136 = cos(t138);
	t135 = sin(t138);
	t133 = t140 * t150 + t149;
	t132 = -t140 * t148 + t151;
	t131 = t134 * t141 + t137 * t144 + pkin(1);
	t130 = t139 * t147 + t140 * t146;
	t1 = [t133 * t135 + t136 * t154, t133 * t136 - t135 * t154, -t140 * t151 + t148, t130 * t142 + t131 * t145 + 0; t132 * t135 - t136 * t152, t132 * t136 + t135 * t152, t140 * t149 + t150, -t130 * t145 + t131 * t142 + 0; -t135 * t153 + t136 * t140, -t135 * t140 - t136 * t153, t139 * t141, -t139 * t146 + t140 * t147 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:01
	% EndTime: 2020-11-04 22:12:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (112->37), mult. (148->56), div. (0->0), fcn. (182->14), ass. (0->39)
	t170 = sin(pkin(6));
	t175 = sin(qJ(2));
	t192 = t170 * t175;
	t176 = sin(qJ(1));
	t191 = t170 * t176;
	t179 = cos(qJ(2));
	t190 = t170 * t179;
	t180 = cos(qJ(1));
	t189 = t170 * t180;
	t188 = t176 * t175;
	t187 = t176 * t179;
	t186 = t180 * t175;
	t185 = t180 * t179;
	t169 = sin(pkin(11));
	t171 = cos(pkin(11));
	t163 = t171 * pkin(5) + t169 * pkin(10) + pkin(4);
	t164 = -t169 * pkin(5) + t171 * pkin(10);
	t174 = sin(qJ(4));
	t178 = cos(qJ(4));
	t157 = t163 * t174 - t164 * t178 + qJ(3);
	t167 = qJ(5) + pkin(2) + pkin(9);
	t184 = t157 * t179 - t167 * t175;
	t172 = cos(pkin(6));
	t159 = t172 * t185 - t188;
	t168 = qJ(4) + pkin(11);
	t165 = sin(t168);
	t166 = cos(t168);
	t183 = t159 * t165 + t166 * t189;
	t161 = t172 * t187 + t186;
	t182 = t161 * t165 + t166 * t191;
	t181 = t163 * t178 + t164 * t174 + pkin(3) + pkin(8);
	t177 = cos(qJ(6));
	t173 = sin(qJ(6));
	t162 = -t172 * t188 + t185;
	t160 = t172 * t186 + t187;
	t158 = -t165 * t190 + t172 * t166;
	t156 = t157 * t175 + t167 * t179 + pkin(1);
	t155 = t181 * t170 + t184 * t172;
	t1 = [t162 * t173 + t182 * t177, t162 * t177 - t182 * t173, -t161 * t166 + t165 * t191, t155 * t176 + t156 * t180 + 0; t160 * t173 - t183 * t177, t160 * t177 + t183 * t173, t159 * t166 - t165 * t189, -t155 * t180 + t156 * t176 + 0; t158 * t177 + t173 * t192, -t158 * t173 + t177 * t192, t172 * t165 + t166 * t190, -t184 * t170 + t181 * t172 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end