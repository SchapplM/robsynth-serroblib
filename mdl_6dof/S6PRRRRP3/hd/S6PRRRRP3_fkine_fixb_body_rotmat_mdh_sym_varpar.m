% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:39
	% EndTime: 2020-11-04 21:18:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:39
	% EndTime: 2020-11-04 21:18:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(pkin(11));
	t98 = sin(pkin(11));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:39
	% EndTime: 2020-11-04 21:18:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t100 = sin(pkin(11));
	t101 = sin(pkin(6));
	t109 = t100 * t101;
	t102 = cos(pkin(11));
	t108 = t102 * t101;
	t103 = cos(pkin(6));
	t104 = sin(qJ(2));
	t107 = t103 * t104;
	t105 = cos(qJ(2));
	t106 = t103 * t105;
	t1 = [-t100 * t107 + t102 * t105, -t100 * t106 - t102 * t104, t109, t102 * pkin(1) + pkin(7) * t109 + 0; t100 * t105 + t102 * t107, -t100 * t104 + t102 * t106, -t108, t100 * pkin(1) - pkin(7) * t108 + 0; t101 * t104, t101 * t105, t103, t103 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:39
	% EndTime: 2020-11-04 21:18:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t113 = sin(pkin(6));
	t126 = t113 * pkin(7);
	t112 = sin(pkin(11));
	t115 = cos(pkin(6));
	t125 = t112 * t115;
	t116 = sin(qJ(3));
	t124 = t113 * t116;
	t118 = cos(qJ(3));
	t123 = t113 * t118;
	t114 = cos(pkin(11));
	t122 = t114 * t115;
	t117 = sin(qJ(2));
	t121 = t115 * t117;
	t119 = cos(qJ(2));
	t120 = t115 * t119;
	t111 = t112 * t119 + t114 * t121;
	t110 = t112 * t121 - t114 * t119;
	t1 = [-t110 * t118 + t112 * t124, t110 * t116 + t112 * t123, t112 * t120 + t114 * t117, (t114 * pkin(2) + pkin(8) * t125) * t119 + (-pkin(2) * t125 + t114 * pkin(8)) * t117 + t112 * t126 + t114 * pkin(1) + 0; t111 * t118 - t114 * t124, -t111 * t116 - t114 * t123, t112 * t117 - t114 * t120, (t112 * pkin(2) - pkin(8) * t122) * t119 + (pkin(2) * t122 + t112 * pkin(8)) * t117 - t114 * t126 + t112 * pkin(1) + 0; t115 * t116 + t117 * t123, t115 * t118 - t117 * t124, -t113 * t119, t115 * pkin(7) + qJ(1) + 0 + (pkin(2) * t117 - pkin(8) * t119) * t113; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:39
	% EndTime: 2020-11-04 21:18:40
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t134 = sin(pkin(6));
	t136 = cos(pkin(6));
	t139 = sin(qJ(2));
	t142 = cos(qJ(2));
	t138 = sin(qJ(3));
	t141 = cos(qJ(3));
	t146 = t141 * pkin(3) + t138 * pkin(9) + pkin(2);
	t144 = -pkin(8) * t142 + t146 * t139;
	t145 = t138 * pkin(3) - t141 * pkin(9) + pkin(7);
	t157 = t145 * t134 - t144 * t136;
	t153 = t134 * t141;
	t152 = t134 * t142;
	t151 = t136 * t139;
	t150 = t136 * t142;
	t149 = t138 * t134;
	t148 = t139 * t141;
	t147 = t141 * t142;
	t143 = pkin(8) * t139 + t146 * t142 + pkin(1);
	t140 = cos(qJ(4));
	t137 = sin(qJ(4));
	t135 = cos(pkin(11));
	t133 = sin(pkin(11));
	t132 = t136 * t148 - t149;
	t131 = t134 * t148 + t136 * t138;
	t130 = t133 * t150 + t135 * t139;
	t129 = t133 * t139 - t135 * t150;
	t128 = -t133 * t132 + t135 * t147;
	t127 = t135 * t132 + t133 * t147;
	t1 = [t128 * t140 + t130 * t137, -t128 * t137 + t130 * t140, -(t133 * t151 - t135 * t142) * t138 - t133 * t153, t157 * t133 + t143 * t135 + 0; t127 * t140 + t129 * t137, -t127 * t137 + t129 * t140, (t133 * t142 + t135 * t151) * t138 + t135 * t153, t143 * t133 - t157 * t135 + 0; t131 * t140 - t137 * t152, -t131 * t137 - t140 * t152, -t136 * t141 + t139 * t149, t144 * t134 + t145 * t136 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:40
	% EndTime: 2020-11-04 21:18:40
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (80->45), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->33)
	t190 = sin(qJ(4)) * pkin(4);
	t173 = sin(pkin(6));
	t189 = t173 * pkin(7);
	t172 = sin(pkin(11));
	t175 = cos(pkin(6));
	t188 = t172 * t175;
	t177 = sin(qJ(3));
	t187 = t173 * t177;
	t179 = cos(qJ(3));
	t186 = t173 * t179;
	t180 = cos(qJ(2));
	t185 = t173 * t180;
	t174 = cos(pkin(11));
	t184 = t174 * t175;
	t178 = sin(qJ(2));
	t183 = t175 * t178;
	t182 = t175 * t180;
	t181 = -pkin(10) - pkin(9);
	t171 = qJ(4) + qJ(5);
	t170 = cos(t171);
	t169 = sin(t171);
	t168 = cos(qJ(4)) * pkin(4) + pkin(3);
	t167 = t175 * t177 + t178 * t186;
	t166 = -t175 * t179 + t178 * t187;
	t165 = t172 * t182 + t174 * t178;
	t164 = t172 * t180 + t174 * t183;
	t163 = t172 * t178 - t174 * t182;
	t162 = t172 * t183 - t174 * t180;
	t161 = -t162 * t179 + t172 * t187;
	t160 = t164 * t179 - t174 * t187;
	t159 = t164 * t177 + t174 * t186;
	t158 = t162 * t177 + t172 * t186;
	t1 = [t161 * t170 + t165 * t169, -t161 * t169 + t165 * t170, -t158, t161 * t168 + t158 * t181 + t165 * t190 + (t174 * pkin(2) + pkin(8) * t188) * t180 + (-pkin(2) * t188 + t174 * pkin(8)) * t178 + t172 * t189 + t174 * pkin(1) + 0; t160 * t170 + t163 * t169, -t160 * t169 + t163 * t170, t159, t160 * t168 - t159 * t181 + t163 * t190 + (t172 * pkin(2) - pkin(8) * t184) * t180 + (pkin(2) * t184 + t172 * pkin(8)) * t178 - t174 * t189 + t172 * pkin(1) + 0; t167 * t170 - t169 * t185, -t167 * t169 - t170 * t185, t166, t175 * pkin(7) - t166 * t181 + t167 * t168 + qJ(1) + 0 + (pkin(2) * t178 + (-pkin(8) - t190) * t180) * t173; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:18:40
	% EndTime: 2020-11-04 21:18:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (95->47), mult. (156->72), div. (0->0), fcn. (202->12), ass. (0->33)
	t208 = sin(pkin(6));
	t222 = t208 * pkin(7);
	t207 = sin(pkin(11));
	t210 = cos(pkin(6));
	t221 = t207 * t210;
	t211 = sin(qJ(3));
	t220 = t208 * t211;
	t213 = cos(qJ(3));
	t219 = t208 * t213;
	t214 = cos(qJ(2));
	t218 = t208 * t214;
	t209 = cos(pkin(11));
	t217 = t209 * t210;
	t212 = sin(qJ(2));
	t216 = t210 * t212;
	t215 = t210 * t214;
	t206 = qJ(4) + qJ(5);
	t205 = -qJ(6) - pkin(10) - pkin(9);
	t204 = cos(t206);
	t203 = sin(t206);
	t202 = pkin(5) * t203 + sin(qJ(4)) * pkin(4);
	t201 = pkin(5) * t204 + cos(qJ(4)) * pkin(4) + pkin(3);
	t200 = t210 * t211 + t212 * t219;
	t199 = -t210 * t213 + t212 * t220;
	t198 = t207 * t215 + t209 * t212;
	t197 = t207 * t214 + t209 * t216;
	t196 = t207 * t212 - t209 * t215;
	t195 = t207 * t216 - t209 * t214;
	t194 = -t195 * t213 + t207 * t220;
	t193 = t197 * t213 - t209 * t220;
	t192 = t197 * t211 + t209 * t219;
	t191 = t195 * t211 + t207 * t219;
	t1 = [t194 * t204 + t198 * t203, -t194 * t203 + t198 * t204, -t191, t194 * t201 + t191 * t205 + t198 * t202 + (t209 * pkin(2) + pkin(8) * t221) * t214 + (-pkin(2) * t221 + t209 * pkin(8)) * t212 + t207 * t222 + t209 * pkin(1) + 0; t193 * t204 + t196 * t203, -t193 * t203 + t196 * t204, t192, t193 * t201 - t192 * t205 + t196 * t202 + (t207 * pkin(2) - pkin(8) * t217) * t214 + (pkin(2) * t217 + t207 * pkin(8)) * t212 - t209 * t222 + t207 * pkin(1) + 0; t200 * t204 - t203 * t218, -t200 * t203 - t204 * t218, t199, t210 * pkin(7) - t199 * t205 + t200 * t201 + qJ(1) + 0 + (pkin(2) * t212 + (-pkin(8) - t202) * t214) * t208; 0, 0, 0, 1;];
	Tc_mdh = t1;
end