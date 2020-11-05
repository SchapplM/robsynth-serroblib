% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:59
	% EndTime: 2020-11-04 22:21:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:59
	% EndTime: 2020-11-04 22:21:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t92 = cos(qJ(1));
	t91 = sin(qJ(1));
	t1 = [t92, -t91, 0, 0; t91, t92, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:59
	% EndTime: 2020-11-04 22:21:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t96 = cos(qJ(1));
	t95 = cos(qJ(2));
	t94 = sin(qJ(1));
	t93 = sin(qJ(2));
	t1 = [t96 * t95, -t96 * t93, t94, t96 * pkin(1) + t94 * pkin(8) + 0; t94 * t95, -t94 * t93, -t96, t94 * pkin(1) - t96 * pkin(8) + 0; t93, t95, 0, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:59
	% EndTime: 2020-11-04 22:21:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t103 = cos(qJ(1));
	t98 = sin(qJ(3));
	t106 = t103 * t98;
	t100 = sin(qJ(1));
	t102 = cos(qJ(2));
	t105 = t100 * t102;
	t101 = cos(qJ(3));
	t104 = t103 * t101;
	t99 = sin(qJ(2));
	t97 = t102 * pkin(2) + t99 * pkin(9) + pkin(1);
	t1 = [t100 * t98 + t102 * t104, t100 * t101 - t102 * t106, t103 * t99, t100 * pkin(8) + t97 * t103 + 0; t101 * t105 - t106, -t98 * t105 - t104, t100 * t99, -t103 * pkin(8) + t97 * t100 + 0; t99 * t101, -t99 * t98, -t102, t99 * pkin(2) - t102 * pkin(9) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:59
	% EndTime: 2020-11-04 22:22:00
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->32), mult. (101->55), div. (0->0), fcn. (135->10), ass. (0->28)
	t116 = sin(pkin(6));
	t133 = qJ(4) * t116;
	t115 = sin(pkin(10));
	t119 = sin(qJ(3));
	t132 = t115 * t119;
	t118 = cos(pkin(6));
	t122 = cos(qJ(3));
	t131 = t118 * t122;
	t117 = cos(pkin(10));
	t130 = t119 * t117;
	t123 = cos(qJ(2));
	t129 = t119 * t123;
	t120 = sin(qJ(2));
	t128 = t120 * t116;
	t127 = t120 * t118;
	t126 = t120 * t122;
	t125 = -t116 * t123 - t119 * t127;
	t124 = cos(qJ(1));
	t121 = sin(qJ(1));
	t114 = qJ(4) * t118 + pkin(9);
	t113 = -pkin(3) * t119 + t122 * t133 - pkin(8);
	t112 = pkin(3) * t122 + t119 * t133 + pkin(2);
	t111 = -t117 * t131 + t132;
	t110 = t115 * t131 + t130;
	t109 = t112 * t123 + t114 * t120 + pkin(1);
	t108 = (-t115 * t122 - t118 * t130) * t123 + t117 * t128;
	t107 = (t117 * t122 - t118 * t132) * t123 + t115 * t128;
	t1 = [t107 * t124 + t110 * t121, t108 * t124 - t111 * t121, (-t121 * t122 + t124 * t129) * t116 + t124 * t127, t109 * t124 - t113 * t121 + 0; t107 * t121 - t110 * t124, t108 * t121 + t111 * t124, (t121 * t129 + t122 * t124) * t116 + t121 * t127, t109 * t121 + t113 * t124 + 0; t115 * t125 + t117 * t126, -t115 * t126 + t117 * t125, -t118 * t123 + t119 * t128, t112 * t120 - t114 * t123 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:00
	% EndTime: 2020-11-04 22:22:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (75->37), mult. (138->61), div. (0->0), fcn. (172->10), ass. (0->30)
	t144 = sin(pkin(10));
	t146 = cos(pkin(10));
	t143 = -t144 * pkin(4) + qJ(5) * t146;
	t145 = sin(pkin(6));
	t147 = cos(pkin(6));
	t139 = t145 * qJ(4) + t143 * t147;
	t142 = pkin(4) * t146 + qJ(5) * t144 + pkin(3);
	t148 = sin(qJ(3));
	t151 = cos(qJ(3));
	t164 = -t139 * t151 + t142 * t148 + pkin(8);
	t162 = t144 * t148;
	t161 = t147 * t151;
	t160 = t148 * t146;
	t152 = cos(qJ(2));
	t159 = t148 * t152;
	t149 = sin(qJ(2));
	t158 = t149 * t145;
	t157 = t149 * t147;
	t156 = t149 * t151;
	t154 = t145 * t152 + t148 * t157;
	t153 = cos(qJ(1));
	t150 = sin(qJ(1));
	t141 = -t146 * t161 + t162;
	t140 = t144 * t161 + t160;
	t138 = -t147 * qJ(4) + t143 * t145 - pkin(9);
	t137 = (t144 * t151 + t147 * t160) * t152 - t146 * t158;
	t136 = (-t146 * t151 + t147 * t162) * t152 - t144 * t158;
	t135 = t139 * t148 + t142 * t151 + pkin(2);
	t134 = t135 * t152 - t138 * t149 + pkin(1);
	t1 = [(-t150 * t151 + t153 * t159) * t145 + t153 * t157, t136 * t153 - t150 * t140, t137 * t153 + t150 * t141, t134 * t153 + t164 * t150 + 0; (t150 * t159 + t153 * t151) * t145 + t150 * t157, t136 * t150 + t153 * t140, t137 * t150 - t153 * t141, t134 * t150 - t164 * t153 + 0; -t152 * t147 + t148 * t158, t154 * t144 - t146 * t156, t144 * t156 + t154 * t146, t135 * t149 + t138 * t152 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:22:00
	% EndTime: 2020-11-04 22:22:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (96->39), mult. (137->61), div. (0->0), fcn. (171->10), ass. (0->32)
	t174 = sin(pkin(10));
	t176 = cos(pkin(10));
	t178 = qJ(6) + pkin(4);
	t173 = qJ(5) * t176 - t178 * t174;
	t175 = sin(pkin(6));
	t177 = cos(pkin(6));
	t179 = qJ(4) + pkin(5);
	t169 = t173 * t177 + t175 * t179;
	t172 = qJ(5) * t174 + t178 * t176 + pkin(3);
	t180 = sin(qJ(3));
	t183 = cos(qJ(3));
	t202 = t169 * t183 - t172 * t180 - pkin(8);
	t199 = -t173 * t175 + t179 * t177 + pkin(9);
	t195 = t174 * t180;
	t193 = t177 * t183;
	t192 = t180 * t176;
	t184 = cos(qJ(2));
	t191 = t180 * t184;
	t181 = sin(qJ(2));
	t190 = t181 * t175;
	t189 = t181 * t177;
	t188 = t181 * t183;
	t186 = t175 * t184 + t180 * t189;
	t185 = cos(qJ(1));
	t182 = sin(qJ(1));
	t171 = -t176 * t193 + t195;
	t170 = t174 * t193 + t192;
	t168 = (t174 * t183 + t177 * t192) * t184 - t176 * t190;
	t167 = (t176 * t183 - t177 * t195) * t184 + t174 * t190;
	t166 = t169 * t180 + t172 * t183 + pkin(2);
	t165 = t166 * t184 + t199 * t181 + pkin(1);
	t1 = [(-t182 * t183 + t185 * t191) * t175 + t185 * t189, t168 * t185 + t182 * t171, t167 * t185 + t182 * t170, t165 * t185 - t202 * t182 + 0; (t182 * t191 + t185 * t183) * t175 + t182 * t189, t168 * t182 - t185 * t171, t167 * t182 - t185 * t170, t165 * t182 + t202 * t185 + 0; -t184 * t177 + t180 * t190, t174 * t188 + t186 * t176, -t186 * t174 + t176 * t188, t166 * t181 - t199 * t184 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end