% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t108 = cos(pkin(10));
	t107 = sin(pkin(10));
	t1 = [t108, -t107, 0, 0; t107, t108, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t109 = sin(pkin(10));
	t110 = sin(pkin(6));
	t118 = t109 * t110;
	t111 = cos(pkin(10));
	t117 = t111 * t110;
	t112 = cos(pkin(6));
	t113 = sin(qJ(2));
	t116 = t112 * t113;
	t114 = cos(qJ(2));
	t115 = t112 * t114;
	t1 = [-t109 * t116 + t111 * t114, -t109 * t115 - t111 * t113, t118, t111 * pkin(1) + pkin(7) * t118 + 0; t109 * t114 + t111 * t116, -t109 * t113 + t111 * t115, -t117, t109 * pkin(1) - pkin(7) * t117 + 0; t110 * t113, t110 * t114, t112, t112 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t122 = sin(pkin(6));
	t135 = t122 * pkin(7);
	t121 = sin(pkin(10));
	t124 = cos(pkin(6));
	t134 = t121 * t124;
	t125 = sin(qJ(3));
	t133 = t122 * t125;
	t127 = cos(qJ(3));
	t132 = t122 * t127;
	t123 = cos(pkin(10));
	t131 = t123 * t124;
	t126 = sin(qJ(2));
	t130 = t124 * t126;
	t128 = cos(qJ(2));
	t129 = t124 * t128;
	t120 = t121 * t128 + t123 * t130;
	t119 = t121 * t130 - t123 * t128;
	t1 = [-t119 * t127 + t121 * t133, t119 * t125 + t121 * t132, t121 * t129 + t123 * t126, (t123 * pkin(2) + pkin(8) * t134) * t128 + (-pkin(2) * t134 + t123 * pkin(8)) * t126 + t121 * t135 + t123 * pkin(1) + 0; t120 * t127 - t123 * t133, -t120 * t125 - t123 * t132, t121 * t126 - t123 * t129, (t121 * pkin(2) - pkin(8) * t131) * t128 + (pkin(2) * t131 + t121 * pkin(8)) * t126 - t123 * t135 + t121 * pkin(1) + 0; t124 * t125 + t126 * t132, t124 * t127 - t126 * t133, -t122 * t128, t124 * pkin(7) + qJ(1) + 0 + (pkin(2) * t126 - pkin(8) * t128) * t122; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t139 = sin(pkin(6));
	t141 = cos(pkin(6));
	t143 = sin(qJ(2));
	t145 = cos(qJ(2));
	t142 = sin(qJ(3));
	t144 = cos(qJ(3));
	t149 = pkin(3) * t144 + qJ(4) * t142 + pkin(2);
	t147 = -pkin(8) * t145 + t149 * t143;
	t148 = t142 * pkin(3) - qJ(4) * t144 + pkin(7);
	t157 = t148 * t139 - t147 * t141;
	t153 = t141 * t143;
	t152 = t141 * t145;
	t151 = t142 * t139;
	t150 = t144 * t139;
	t146 = pkin(8) * t143 + t149 * t145 + pkin(1);
	t140 = cos(pkin(10));
	t138 = sin(pkin(10));
	t137 = t138 * t145 + t140 * t153;
	t136 = t138 * t153 - t140 * t145;
	t1 = [-t136 * t144 + t138 * t151, t138 * t152 + t140 * t143, -t136 * t142 - t138 * t150, t157 * t138 + t146 * t140 + 0; t137 * t144 - t140 * t151, t138 * t143 - t140 * t152, t137 * t142 + t140 * t150, t146 * t138 - t157 * t140 + 0; t141 * t142 + t143 * t150, -t139 * t145, -t141 * t144 + t143 * t151, t147 * t139 + t148 * t141 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (57->27), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->22)
	t161 = sin(pkin(6));
	t163 = cos(pkin(6));
	t164 = qJ(5) - pkin(8);
	t166 = sin(qJ(2));
	t168 = cos(qJ(2));
	t165 = sin(qJ(3));
	t167 = cos(qJ(3));
	t169 = pkin(3) + pkin(4);
	t173 = qJ(4) * t165 + t169 * t167 + pkin(2);
	t171 = t164 * t168 + t173 * t166;
	t172 = qJ(4) * t167 - t169 * t165 - pkin(7);
	t181 = t172 * t161 + t171 * t163;
	t178 = t161 * t165;
	t177 = t163 * t166;
	t176 = t163 * t168;
	t174 = t167 * t161;
	t170 = -t164 * t166 + t173 * t168 + pkin(1);
	t162 = cos(pkin(10));
	t160 = sin(pkin(10));
	t159 = t160 * t168 + t162 * t177;
	t158 = t160 * t177 - t162 * t168;
	t1 = [-t158 * t165 - t160 * t174, t158 * t167 - t160 * t178, -t160 * t176 - t162 * t166, -t181 * t160 + t170 * t162 + 0; t159 * t165 + t162 * t174, -t159 * t167 + t162 * t178, -t160 * t166 + t162 * t176, t170 * t160 + t181 * t162 + 0; -t163 * t167 + t166 * t178, -t163 * t165 - t166 * t174, t161 * t168, t171 * t161 - t172 * t163 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:15
	% EndTime: 2020-11-04 21:07:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (86->38), mult. (134->60), div. (0->0), fcn. (168->10), ass. (0->33)
	t189 = pkin(3) + pkin(4) + pkin(9);
	t191 = sin(pkin(6));
	t193 = cos(pkin(6));
	t195 = qJ(4) + pkin(5);
	t198 = sin(qJ(2));
	t200 = cos(qJ(3));
	t197 = sin(qJ(3));
	t204 = t189 * t197 + pkin(7);
	t205 = t195 * t197 + pkin(2);
	t194 = qJ(5) - pkin(8);
	t201 = cos(qJ(2));
	t209 = t194 * t201;
	t211 = t193 * t198;
	t217 = (t205 * t198 + t209) * t193 - t204 * t191 + (t189 * t211 + t191 * t195) * t200;
	t213 = t191 * t197;
	t212 = t191 * t201;
	t210 = t193 * t201;
	t208 = t197 * t198;
	t207 = t197 * t201;
	t206 = t200 * t191;
	t203 = t189 * t200 + t205;
	t202 = -t194 * t198 + t203 * t201 + pkin(1);
	t199 = cos(qJ(6));
	t196 = sin(qJ(6));
	t192 = cos(pkin(10));
	t190 = sin(pkin(10));
	t188 = t191 * t208 - t193 * t200;
	t187 = t193 * t208 + t206;
	t186 = t190 * t210 + t192 * t198;
	t185 = t190 * t198 - t192 * t210;
	t183 = -t190 * t187 + t192 * t207;
	t182 = t192 * t187 + t190 * t207;
	t1 = [t183 * t199 - t196 * t186, -t183 * t196 - t199 * t186, -(t190 * t211 - t192 * t201) * t200 + t190 * t213, -t217 * t190 + t202 * t192 + 0; t182 * t199 - t196 * t185, -t182 * t196 - t199 * t185, (t190 * t201 + t192 * t211) * t200 - t192 * t213, t202 * t190 + t217 * t192 + 0; t188 * t199 + t196 * t212, -t188 * t196 + t199 * t212, t193 * t197 + t198 * t206, qJ(1) + 0 + (t203 * t198 + t209) * t191 + (-t195 * t200 + t204) * t193; 0, 0, 0, 1;];
	Tc_mdh = t1;
end