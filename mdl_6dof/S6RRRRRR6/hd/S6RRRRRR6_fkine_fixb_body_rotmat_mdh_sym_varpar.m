% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t82 = cos(qJ(1));
	t81 = sin(qJ(1));
	t1 = [t82, -t81, 0, 0; t81, t82, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t83 = sin(pkin(6));
	t86 = sin(qJ(1));
	t94 = t86 * t83;
	t85 = sin(qJ(2));
	t93 = t86 * t85;
	t87 = cos(qJ(2));
	t92 = t86 * t87;
	t88 = cos(qJ(1));
	t91 = t88 * t83;
	t90 = t88 * t85;
	t89 = t88 * t87;
	t84 = cos(pkin(6));
	t1 = [-t84 * t93 + t89, -t84 * t92 - t90, t94, t88 * pkin(1) + pkin(8) * t94 + 0; t84 * t90 + t92, t84 * t89 - t93, -t91, t86 * pkin(1) - pkin(8) * t91 + 0; t83 * t85, t83 * t87, t84, t84 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t100 = sin(qJ(3));
	t98 = sin(pkin(6));
	t113 = t100 * t98;
	t103 = cos(qJ(3));
	t112 = t98 * t103;
	t101 = sin(qJ(2));
	t111 = t101 * t103;
	t102 = sin(qJ(1));
	t110 = t102 * t101;
	t104 = cos(qJ(2));
	t109 = t102 * t104;
	t105 = cos(qJ(1));
	t108 = t105 * t101;
	t107 = t105 * t104;
	t106 = pkin(2) * t101 - pkin(9) * t104;
	t99 = cos(pkin(6));
	t97 = t104 * pkin(2) + t101 * pkin(9) + pkin(1);
	t96 = -t99 * t108 - t109;
	t95 = t98 * pkin(8) - t106 * t99;
	t1 = [(-t99 * t111 + t113) * t102 + t103 * t107, (t99 * t110 - t107) * t100 + t102 * t112, t99 * t109 + t108, t95 * t102 + t97 * t105 + 0; -t96 * t103 - t105 * t113, t96 * t100 - t105 * t112, -t99 * t107 + t110, t97 * t102 - t95 * t105 + 0; t99 * t100 + t98 * t111, -t101 * t113 + t99 * t103, -t98 * t104, t99 * pkin(8) + t106 * t98 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->26), mult. (69->38), div. (0->0), fcn. (90->10), ass. (0->25)
	t122 = sin(pkin(6));
	t125 = sin(qJ(2));
	t138 = t122 * t125;
	t126 = sin(qJ(1));
	t137 = t122 * t126;
	t128 = cos(qJ(1));
	t136 = t122 * t128;
	t135 = t126 * t125;
	t127 = cos(qJ(2));
	t134 = t126 * t127;
	t133 = t128 * t125;
	t132 = t128 * t127;
	t131 = sin(qJ(3)) * pkin(3) + pkin(8);
	t118 = cos(qJ(3)) * pkin(3) + pkin(2);
	t129 = pkin(10) + pkin(9);
	t130 = t118 * t125 - t129 * t127;
	t123 = cos(pkin(6));
	t121 = qJ(3) + qJ(4);
	t120 = cos(t121);
	t119 = sin(t121);
	t117 = t123 * t133 + t134;
	t116 = t123 * t135 - t132;
	t115 = t118 * t127 + t129 * t125 + pkin(1);
	t114 = t122 * t131 - t130 * t123;
	t1 = [-t116 * t120 + t119 * t137, t116 * t119 + t120 * t137, t123 * t134 + t133, t114 * t126 + t115 * t128 + 0; t117 * t120 - t119 * t136, -t117 * t119 - t120 * t136, -t123 * t132 + t135, -t114 * t128 + t115 * t126 + 0; t123 * t119 + t120 * t138, -t119 * t138 + t123 * t120, -t122 * t127, t130 * t122 + t131 * t123 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (100->38), mult. (139->56), div. (0->0), fcn. (183->12), ass. (0->36)
	t155 = sin(pkin(6));
	t159 = sin(qJ(2));
	t174 = t155 * t159;
	t160 = sin(qJ(1));
	t173 = t155 * t160;
	t162 = cos(qJ(2));
	t172 = t155 * t162;
	t163 = cos(qJ(1));
	t171 = t155 * t163;
	t170 = t160 * t159;
	t169 = t160 * t162;
	t168 = t163 * t159;
	t167 = t163 * t162;
	t166 = sin(qJ(3)) * pkin(3) + pkin(8);
	t151 = cos(qJ(3)) * pkin(3) + pkin(2);
	t164 = pkin(10) + pkin(9);
	t165 = t151 * t159 - t164 * t162;
	t161 = cos(qJ(5));
	t157 = sin(qJ(5));
	t156 = cos(pkin(6));
	t154 = qJ(3) + qJ(4);
	t153 = cos(t154);
	t152 = sin(t154);
	t150 = t156 * t169 + t168;
	t149 = t156 * t168 + t169;
	t148 = -t156 * t167 + t170;
	t147 = t156 * t170 - t167;
	t146 = t151 * t162 + t164 * t159 + pkin(1);
	t145 = t156 * t152 + t153 * t174;
	t144 = t152 * t174 - t156 * t153;
	t143 = t155 * t166 - t165 * t156;
	t142 = -t147 * t153 + t152 * t173;
	t141 = t149 * t153 - t152 * t171;
	t140 = t149 * t152 + t153 * t171;
	t139 = t147 * t152 + t153 * t173;
	t1 = [t142 * t161 + t150 * t157, -t142 * t157 + t150 * t161, -t139, t142 * pkin(4) - t139 * pkin(11) + t143 * t160 + t146 * t163 + 0; t141 * t161 + t148 * t157, -t141 * t157 + t148 * t161, t140, t141 * pkin(4) + t140 * pkin(11) - t143 * t163 + t146 * t160 + 0; t145 * t161 - t157 * t172, -t145 * t157 - t161 * t172, t144, t145 * pkin(4) + t144 * pkin(11) + t165 * t155 + t166 * t156 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:47:44
	% EndTime: 2020-11-04 22:47:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (123->45), mult. (155->61), div. (0->0), fcn. (201->14), ass. (0->40)
	t215 = pkin(5) * sin(qJ(5));
	t188 = cos(qJ(3)) * pkin(3) + pkin(2);
	t199 = sin(qJ(2));
	t214 = t188 * t199;
	t195 = sin(pkin(6));
	t213 = t195 * t199;
	t200 = sin(qJ(1));
	t212 = t195 * t200;
	t201 = cos(qJ(2));
	t211 = t195 * t201;
	t202 = cos(qJ(1));
	t210 = t195 * t202;
	t209 = t200 * t199;
	t208 = t200 * t201;
	t207 = t202 * t199;
	t206 = t202 * t201;
	t205 = sin(qJ(3)) * pkin(3) + pkin(8);
	t204 = pkin(10) + pkin(9);
	t203 = -pkin(12) - pkin(11);
	t196 = cos(pkin(6));
	t194 = qJ(3) + qJ(4);
	t193 = qJ(5) + qJ(6);
	t192 = cos(t194);
	t191 = cos(t193);
	t190 = sin(t194);
	t189 = sin(t193);
	t187 = cos(qJ(5)) * pkin(5) + pkin(4);
	t186 = t196 * t208 + t207;
	t185 = t196 * t207 + t208;
	t184 = -t196 * t206 + t209;
	t183 = t196 * t209 - t206;
	t182 = t188 * t201 + t204 * t199 + pkin(1);
	t181 = t196 * t190 + t192 * t213;
	t180 = t190 * t213 - t196 * t192;
	t179 = t195 * t205 + (t204 * t201 - t214) * t196;
	t178 = -t183 * t192 + t190 * t212;
	t177 = t185 * t192 - t190 * t210;
	t176 = t185 * t190 + t192 * t210;
	t175 = t183 * t190 + t192 * t212;
	t1 = [t178 * t191 + t186 * t189, -t178 * t189 + t186 * t191, -t175, t175 * t203 + t178 * t187 + t179 * t200 + t182 * t202 + t186 * t215 + 0; t177 * t191 + t184 * t189, -t177 * t189 + t184 * t191, t176, -t176 * t203 + t177 * t187 - t179 * t202 + t182 * t200 + t184 * t215 + 0; t181 * t191 - t189 * t211, -t181 * t189 - t191 * t211, t180, -t180 * t203 + t181 * t187 + pkin(7) + 0 + t205 * t196 + (t214 + (-t204 - t215) * t201) * t195; 0, 0, 0, 1;];
	Tc_mdh = t1;
end