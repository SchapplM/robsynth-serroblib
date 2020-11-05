% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t1 = [t83, -t82, 0, 0; t82, t83, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.02s
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
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t101 = sin(qJ(3));
	t99 = sin(pkin(6));
	t114 = t101 * t99;
	t104 = cos(qJ(3));
	t113 = t99 * t104;
	t102 = sin(qJ(2));
	t112 = t102 * t104;
	t103 = sin(qJ(1));
	t111 = t103 * t102;
	t105 = cos(qJ(2));
	t110 = t103 * t105;
	t106 = cos(qJ(1));
	t109 = t106 * t102;
	t108 = t106 * t105;
	t107 = pkin(2) * t102 - pkin(9) * t105;
	t100 = cos(pkin(6));
	t98 = pkin(2) * t105 + pkin(9) * t102 + pkin(1);
	t97 = -t100 * t109 - t110;
	t96 = t99 * pkin(8) - t100 * t107;
	t1 = [(-t100 * t112 + t114) * t103 + t104 * t108, (t100 * t111 - t108) * t101 + t103 * t113, t100 * t110 + t109, t103 * t96 + t106 * t98 + 0; -t104 * t97 - t106 * t114, t101 * t97 - t106 * t113, -t100 * t108 + t111, t103 * t98 - t106 * t96 + 0; t100 * t101 + t112 * t99, t100 * t104 - t102 * t114, -t99 * t105, t100 * pkin(8) + t107 * t99 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t123 = sin(pkin(6));
	t126 = sin(qJ(2));
	t139 = t123 * t126;
	t127 = sin(qJ(1));
	t138 = t123 * t127;
	t129 = cos(qJ(1));
	t137 = t123 * t129;
	t136 = t127 * t126;
	t128 = cos(qJ(2));
	t135 = t127 * t128;
	t134 = t129 * t126;
	t133 = t129 * t128;
	t132 = sin(qJ(3)) * pkin(3) + pkin(8);
	t119 = cos(qJ(3)) * pkin(3) + pkin(2);
	t130 = pkin(10) + pkin(9);
	t131 = t119 * t126 - t130 * t128;
	t124 = cos(pkin(6));
	t122 = qJ(3) + qJ(4);
	t121 = cos(t122);
	t120 = sin(t122);
	t118 = t124 * t134 + t135;
	t117 = t124 * t136 - t133;
	t116 = t119 * t128 + t130 * t126 + pkin(1);
	t115 = t123 * t132 - t131 * t124;
	t1 = [-t117 * t121 + t120 * t138, t117 * t120 + t121 * t138, t124 * t135 + t134, t115 * t127 + t116 * t129 + 0; t118 * t121 - t120 * t137, -t118 * t120 - t121 * t137, -t124 * t133 + t136, -t115 * t129 + t116 * t127 + 0; t124 * t120 + t121 * t139, -t120 * t139 + t124 * t121, -t123 * t128, t131 * t123 + t132 * t124 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (139->56), div. (0->0), fcn. (183->12), ass. (0->36)
	t156 = sin(pkin(6));
	t160 = sin(qJ(2));
	t175 = t156 * t160;
	t161 = sin(qJ(1));
	t174 = t156 * t161;
	t163 = cos(qJ(2));
	t173 = t156 * t163;
	t164 = cos(qJ(1));
	t172 = t156 * t164;
	t171 = t161 * t160;
	t170 = t161 * t163;
	t169 = t164 * t160;
	t168 = t164 * t163;
	t167 = sin(qJ(3)) * pkin(3) + pkin(8);
	t152 = cos(qJ(3)) * pkin(3) + pkin(2);
	t165 = pkin(10) + pkin(9);
	t166 = t152 * t160 - t165 * t163;
	t162 = cos(qJ(5));
	t158 = sin(qJ(5));
	t157 = cos(pkin(6));
	t155 = qJ(3) + qJ(4);
	t154 = cos(t155);
	t153 = sin(t155);
	t151 = t157 * t170 + t169;
	t150 = t157 * t169 + t170;
	t149 = -t157 * t168 + t171;
	t148 = t157 * t171 - t168;
	t147 = t152 * t163 + t165 * t160 + pkin(1);
	t146 = t157 * t153 + t154 * t175;
	t145 = t153 * t175 - t157 * t154;
	t144 = t156 * t167 - t166 * t157;
	t143 = -t148 * t154 + t153 * t174;
	t142 = t150 * t154 - t153 * t172;
	t141 = t150 * t153 + t154 * t172;
	t140 = t148 * t153 + t154 * t174;
	t1 = [t143 * t162 + t151 * t158, -t143 * t158 + t151 * t162, -t140, t143 * pkin(4) - t140 * pkin(11) + t144 * t161 + t147 * t164 + 0; t142 * t162 + t149 * t158, -t142 * t158 + t149 * t162, t141, t142 * pkin(4) + t141 * pkin(11) - t144 * t164 + t147 * t161 + 0; t146 * t162 - t158 * t173, -t146 * t158 - t162 * t173, t145, t146 * pkin(4) + t145 * pkin(11) + t166 * t156 + t167 * t157 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:10
	% EndTime: 2020-11-04 22:44:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (111->44), mult. (155->61), div. (0->0), fcn. (201->12), ass. (0->40)
	t194 = cos(pkin(6));
	t201 = cos(qJ(2));
	t202 = cos(qJ(1));
	t205 = t202 * t201;
	t198 = sin(qJ(2));
	t199 = sin(qJ(1));
	t208 = t199 * t198;
	t185 = -t194 * t205 + t208;
	t196 = sin(qJ(5));
	t215 = t185 * t196;
	t206 = t202 * t198;
	t207 = t199 * t201;
	t187 = t194 * t207 + t206;
	t214 = t187 * t196;
	t189 = cos(qJ(3)) * pkin(3) + pkin(2);
	t213 = t189 * t198;
	t193 = sin(pkin(6));
	t212 = t193 * t198;
	t211 = t193 * t199;
	t210 = t193 * t201;
	t209 = t193 * t202;
	t204 = sin(qJ(3)) * pkin(3) + pkin(8);
	t203 = pkin(10) + pkin(9);
	t200 = cos(qJ(5));
	t195 = -qJ(6) - pkin(11);
	t192 = qJ(3) + qJ(4);
	t191 = cos(t192);
	t190 = sin(t192);
	t188 = t200 * pkin(5) + pkin(4);
	t186 = t194 * t206 + t207;
	t184 = t194 * t208 - t205;
	t183 = t189 * t201 + t203 * t198 + pkin(1);
	t182 = t194 * t190 + t191 * t212;
	t181 = t190 * t212 - t194 * t191;
	t180 = t193 * t204 + (t203 * t201 - t213) * t194;
	t179 = -t184 * t191 + t190 * t211;
	t178 = t186 * t191 - t190 * t209;
	t177 = t186 * t190 + t191 * t209;
	t176 = t184 * t190 + t191 * t211;
	t1 = [t179 * t200 + t214, -t179 * t196 + t187 * t200, -t176, pkin(5) * t214 + t176 * t195 + t179 * t188 + t180 * t199 + t183 * t202 + 0; t178 * t200 + t215, -t178 * t196 + t185 * t200, t177, pkin(5) * t215 - t177 * t195 + t178 * t188 - t180 * t202 + t183 * t199 + 0; t182 * t200 - t196 * t210, -t182 * t196 - t200 * t210, t181, -t181 * t195 + t182 * t188 + pkin(7) + 0 + t204 * t194 + (t213 + (-pkin(5) * t196 - t203) * t201) * t193; 0, 0, 0, 1;];
	Tc_mdh = t1;
end