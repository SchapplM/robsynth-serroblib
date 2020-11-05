% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t89 = cos(qJ(1));
	t88 = sin(qJ(1));
	t1 = [t89, -t88, 0, 0; t88, t89, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t90 = sin(pkin(6));
	t93 = sin(qJ(1));
	t101 = t93 * t90;
	t92 = sin(qJ(2));
	t100 = t93 * t92;
	t94 = cos(qJ(2));
	t99 = t93 * t94;
	t95 = cos(qJ(1));
	t98 = t95 * t90;
	t97 = t95 * t92;
	t96 = t95 * t94;
	t91 = cos(pkin(6));
	t1 = [-t91 * t100 + t96, -t91 * t99 - t97, t101, t95 * pkin(1) + pkin(8) * t101 + 0; t91 * t97 + t99, t91 * t96 - t100, -t98, t93 * pkin(1) - pkin(8) * t98 + 0; t90 * t92, t90 * t94, t91, t91 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t105 = sin(pkin(6));
	t110 = cos(qJ(3));
	t120 = t105 * t110;
	t107 = sin(qJ(3));
	t119 = t107 * t105;
	t108 = sin(qJ(2));
	t118 = t108 * t110;
	t109 = sin(qJ(1));
	t117 = t109 * t108;
	t111 = cos(qJ(2));
	t116 = t109 * t111;
	t112 = cos(qJ(1));
	t115 = t112 * t108;
	t114 = t112 * t111;
	t113 = pkin(2) * t108 - pkin(9) * t111;
	t106 = cos(pkin(6));
	t104 = t111 * pkin(2) + t108 * pkin(9) + pkin(1);
	t103 = t106 * t115 + t116;
	t102 = t105 * pkin(8) - t113 * t106;
	t1 = [(-t106 * t118 + t119) * t109 + t110 * t114, (t106 * t117 - t114) * t107 + t109 * t120, t106 * t116 + t115, t102 * t109 + t104 * t112 + 0; t103 * t110 - t112 * t119, -t103 * t107 - t112 * t120, -t106 * t114 + t117, -t102 * t112 + t104 * t109 + 0; t105 * t118 + t106 * t107, t106 * t110 - t108 * t119, -t105 * t111, t106 * pkin(8) + t113 * t105 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t129 = sin(pkin(6));
	t132 = sin(qJ(2));
	t145 = t129 * t132;
	t133 = sin(qJ(1));
	t144 = t129 * t133;
	t135 = cos(qJ(1));
	t143 = t129 * t135;
	t142 = t133 * t132;
	t134 = cos(qJ(2));
	t141 = t133 * t134;
	t140 = t135 * t132;
	t139 = t135 * t134;
	t138 = sin(qJ(3)) * pkin(3) + pkin(8);
	t125 = cos(qJ(3)) * pkin(3) + pkin(2);
	t136 = pkin(10) + pkin(9);
	t137 = t125 * t132 - t136 * t134;
	t130 = cos(pkin(6));
	t128 = qJ(3) + qJ(4);
	t127 = cos(t128);
	t126 = sin(t128);
	t124 = -t130 * t142 + t139;
	t123 = t130 * t140 + t141;
	t122 = t125 * t134 + t136 * t132 + pkin(1);
	t121 = t129 * t138 - t137 * t130;
	t1 = [t124 * t127 + t126 * t144, -t124 * t126 + t127 * t144, t130 * t141 + t140, t121 * t133 + t122 * t135 + 0; t123 * t127 - t126 * t143, -t123 * t126 - t127 * t143, -t130 * t139 + t142, -t121 * t135 + t122 * t133 + 0; t130 * t126 + t127 * t145, -t126 * t145 + t130 * t127, -t129 * t134, t137 * t129 + t138 * t130 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (82->36), mult. (105->44), div. (0->0), fcn. (136->10), ass. (0->31)
	t160 = sin(pkin(6));
	t163 = sin(qJ(2));
	t176 = t160 * t163;
	t164 = sin(qJ(1));
	t175 = t160 * t164;
	t166 = cos(qJ(1));
	t174 = t160 * t166;
	t173 = t164 * t163;
	t165 = cos(qJ(2));
	t172 = t164 * t165;
	t171 = t166 * t163;
	t170 = t166 * t165;
	t169 = sin(qJ(3)) * pkin(3) + pkin(8);
	t156 = cos(qJ(3)) * pkin(3) + pkin(2);
	t167 = pkin(10) + pkin(9);
	t168 = t156 * t163 - t167 * t165;
	t161 = cos(pkin(6));
	t159 = qJ(3) + qJ(4);
	t158 = cos(t159);
	t157 = sin(t159);
	t155 = -t161 * t173 + t170;
	t154 = t161 * t171 + t172;
	t153 = t156 * t165 + t167 * t163 + pkin(1);
	t152 = t161 * t157 + t158 * t176;
	t151 = t157 * t176 - t161 * t158;
	t150 = t160 * t169 - t168 * t161;
	t149 = -t155 * t157 + t158 * t175;
	t148 = t155 * t158 + t157 * t175;
	t147 = t154 * t158 - t157 * t174;
	t146 = t154 * t157 + t158 * t174;
	t1 = [t161 * t172 + t171, -t148, -t149, t148 * pkin(4) - t149 * qJ(5) + t150 * t164 + t153 * t166 + 0; -t161 * t170 + t173, -t147, t146, t147 * pkin(4) + t146 * qJ(5) - t150 * t166 + t153 * t164 + 0; -t160 * t165, -t152, t151, t152 * pkin(4) + t151 * qJ(5) + t168 * t160 + t169 * t161 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:40:13
	% EndTime: 2020-11-04 22:40:14
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (125->40), mult. (152->61), div. (0->0), fcn. (186->14), ass. (0->41)
	t201 = cos(qJ(4));
	t218 = sin(qJ(4));
	t219 = pkin(11) + pkin(4);
	t187 = qJ(5) * t218 + t219 * t201 + pkin(3);
	t197 = sin(qJ(3));
	t181 = t187 * t197 + pkin(8);
	t188 = qJ(5) * t201 - t218 * t219;
	t182 = -t188 * t197 - pkin(2);
	t194 = sin(pkin(6));
	t195 = cos(pkin(6));
	t198 = sin(qJ(2));
	t202 = cos(qJ(3));
	t192 = pkin(5) + pkin(9) + pkin(10);
	t203 = cos(qJ(2));
	t217 = t192 * t203;
	t220 = (t182 * t198 + t217) * t195 - (t187 * t195 * t198 + t194 * t188) * t202 + t181 * t194;
	t216 = t194 * t198;
	t199 = sin(qJ(1));
	t215 = t194 * t199;
	t214 = t194 * t203;
	t204 = cos(qJ(1));
	t213 = t194 * t204;
	t212 = t199 * t198;
	t211 = t199 * t203;
	t210 = t204 * t198;
	t209 = t204 * t203;
	t184 = t195 * t210 + t211;
	t193 = qJ(3) + qJ(4);
	t190 = sin(t193);
	t191 = cos(t193);
	t206 = t184 * t190 + t191 * t213;
	t186 = -t195 * t212 + t209;
	t205 = t186 * t190 - t191 * t215;
	t200 = cos(qJ(6));
	t196 = sin(qJ(6));
	t185 = t195 * t211 + t210;
	t183 = t195 * t209 - t212;
	t180 = t190 * t216 - t195 * t191;
	t179 = t187 * t202 - t182;
	t177 = t179 * t203 + t192 * t198 + pkin(1);
	t1 = [t185 * t200 + t205 * t196, -t185 * t196 + t205 * t200, t186 * t191 + t190 * t215, t177 * t204 + t220 * t199 + 0; -t200 * t183 + t206 * t196, t196 * t183 + t206 * t200, t184 * t191 - t190 * t213, t177 * t199 - t220 * t204 + 0; t180 * t196 - t200 * t214, t180 * t200 + t196 * t214, t195 * t190 + t191 * t216, pkin(7) + 0 + (t179 * t198 - t217) * t194 + (-t188 * t202 + t181) * t195; 0, 0, 0, 1;];
	Tc_mdh = t1;
end