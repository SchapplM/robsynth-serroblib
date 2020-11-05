% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(pkin(11));
	t98 = sin(pkin(11));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
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
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t129 = pkin(2) * sin(qJ(2));
	t122 = qJ(2) + pkin(12);
	t127 = pkin(7) + qJ(3);
	t126 = cos(pkin(6));
	t125 = cos(pkin(11));
	t124 = sin(pkin(6));
	t123 = sin(pkin(11));
	t121 = cos(t122);
	t120 = sin(t122);
	t119 = pkin(6) - t122;
	t118 = pkin(6) + t122;
	t117 = cos(qJ(2)) * pkin(2) + pkin(1);
	t116 = cos(t118);
	t115 = sin(t119);
	t114 = cos(t119) / 0.2e1;
	t113 = sin(t118) / 0.2e1;
	t112 = -t124 * t127 + t126 * t129;
	t111 = t116 / 0.2e1 + t114;
	t110 = t113 - t115 / 0.2e1;
	t1 = [-t123 * t110 + t125 * t121, -t123 * t111 - t125 * t120, t123 * t124, -t123 * t112 + t125 * t117 + 0; t125 * t110 + t123 * t121, t125 * t111 - t123 * t120, -t125 * t124, t125 * t112 + t123 * t117 + 0; t114 - t116 / 0.2e1, t115 / 0.2e1 + t113, t126, t124 * t129 + t126 * t127 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t141 = sin(pkin(11));
	t145 = cos(pkin(6));
	t155 = t141 * t145;
	t142 = sin(pkin(6));
	t146 = pkin(7) + qJ(3);
	t154 = t142 * t146;
	t147 = sin(qJ(4));
	t153 = t142 * t147;
	t149 = cos(qJ(4));
	t152 = t142 * t149;
	t144 = cos(pkin(11));
	t151 = t144 * t145;
	t139 = qJ(2) + pkin(12);
	t150 = cos(qJ(2));
	t148 = sin(qJ(2));
	t143 = cos(pkin(12));
	t140 = sin(pkin(12));
	t138 = cos(t139);
	t137 = sin(t139);
	t136 = pkin(6) - t139;
	t135 = pkin(6) + t139;
	t134 = -t140 * pkin(3) + t143 * pkin(8);
	t133 = t143 * pkin(3) + t140 * pkin(8) + pkin(2);
	t132 = cos(t135) + cos(t136);
	t131 = t137 * t151 + t141 * t138;
	t130 = t137 * t155 - t144 * t138;
	t1 = [-t130 * t149 + t141 * t153, t130 * t147 + t141 * t152, t144 * t137 + t141 * t132 / 0.2e1, (t144 * t133 + t134 * t155) * t150 + (-t133 * t155 + t144 * t134) * t148 + t141 * t154 + t144 * pkin(1) + 0; t131 * t149 - t144 * t153, -t131 * t147 - t144 * t152, t141 * t137 - t144 * t132 / 0.2e1, (t141 * t133 - t134 * t151) * t150 + (t133 * t151 + t141 * t134) * t148 - t144 * t154 + t141 * pkin(1) + 0; t137 * t152 + t145 * t147, -t137 * t153 + t145 * t149, -sin(t136) / 0.2e1 - sin(t135) / 0.2e1, t145 * t146 + qJ(1) + 0 + (t133 * t148 - t134 * t150) * t142; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (134->48), mult. (156->71), div. (0->0), fcn. (198->18), ass. (0->38)
	t176 = sin(pkin(11));
	t180 = cos(pkin(6));
	t192 = t176 * t180;
	t177 = sin(pkin(6));
	t181 = pkin(7) + qJ(3);
	t191 = t177 * t181;
	t183 = sin(qJ(4));
	t190 = t177 * t183;
	t186 = cos(qJ(4));
	t189 = t177 * t186;
	t179 = cos(pkin(11));
	t188 = t179 * t180;
	t174 = qJ(2) + pkin(12);
	t187 = cos(qJ(2));
	t185 = cos(qJ(5));
	t184 = sin(qJ(2));
	t182 = sin(qJ(5));
	t178 = cos(pkin(12));
	t175 = sin(pkin(12));
	t173 = cos(t174);
	t172 = sin(t174);
	t171 = pkin(6) - t174;
	t170 = pkin(6) + t174;
	t169 = -t175 * pkin(3) + t178 * pkin(8);
	t168 = t178 * pkin(3) + t175 * pkin(8) + pkin(2);
	t167 = cos(t170) + cos(t171);
	t166 = -sin(t171) / 0.2e1 - sin(t170) / 0.2e1;
	t165 = t172 * t189 + t180 * t183;
	t164 = t172 * t190 - t180 * t186;
	t163 = t172 * t188 + t176 * t173;
	t162 = t172 * t192 - t179 * t173;
	t161 = t179 * t172 + t176 * t167 / 0.2e1;
	t160 = t176 * t172 - t179 * t167 / 0.2e1;
	t159 = -t162 * t186 + t176 * t190;
	t158 = t163 * t186 - t179 * t190;
	t157 = t163 * t183 + t179 * t189;
	t156 = t162 * t183 + t176 * t189;
	t1 = [t159 * t185 + t161 * t182, -t159 * t182 + t161 * t185, -t156, t159 * pkin(4) - t156 * pkin(9) + (t179 * t168 + t169 * t192) * t187 + (-t168 * t192 + t179 * t169) * t184 + t176 * t191 + t179 * pkin(1) + 0; t158 * t185 + t160 * t182, -t158 * t182 + t160 * t185, t157, t158 * pkin(4) + t157 * pkin(9) + (t176 * t168 - t169 * t188) * t187 + (t168 * t188 + t176 * t169) * t184 - t179 * t191 + t176 * pkin(1) + 0; t165 * t185 + t166 * t182, -t165 * t182 + t166 * t185, t164, t165 * pkin(4) + t164 * pkin(9) + t180 * t181 + qJ(1) + 0 + (t168 * t184 - t169 * t187) * t177; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:04:04
	% EndTime: 2020-11-04 21:04:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (174->54), mult. (173->76), div. (0->0), fcn. (216->20), ass. (0->42)
	t234 = pkin(5) * sin(qJ(5));
	t217 = sin(pkin(11));
	t221 = cos(pkin(6));
	t233 = t217 * t221;
	t218 = sin(pkin(6));
	t222 = pkin(7) + qJ(3);
	t232 = t218 * t222;
	t224 = sin(qJ(4));
	t231 = t218 * t224;
	t226 = cos(qJ(4));
	t230 = t218 * t226;
	t220 = cos(pkin(11));
	t229 = t220 * t221;
	t214 = qJ(2) + pkin(12);
	t228 = -pkin(10) - pkin(9);
	t227 = cos(qJ(2));
	t225 = sin(qJ(2));
	t219 = cos(pkin(12));
	t216 = sin(pkin(12));
	t215 = qJ(5) + qJ(6);
	t213 = cos(t215);
	t212 = sin(t215);
	t211 = cos(t214);
	t210 = sin(t214);
	t209 = pkin(6) - t214;
	t208 = pkin(6) + t214;
	t207 = cos(qJ(5)) * pkin(5) + pkin(4);
	t206 = -t216 * pkin(3) + t219 * pkin(8);
	t205 = t219 * pkin(3) + t216 * pkin(8) + pkin(2);
	t204 = cos(t208) + cos(t209);
	t203 = -sin(t209) / 0.2e1 - sin(t208) / 0.2e1;
	t202 = t210 * t230 + t221 * t224;
	t201 = t210 * t231 - t221 * t226;
	t200 = t210 * t229 + t217 * t211;
	t199 = t210 * t233 - t220 * t211;
	t198 = t220 * t210 + t217 * t204 / 0.2e1;
	t197 = t217 * t210 - t220 * t204 / 0.2e1;
	t196 = -t199 * t226 + t217 * t231;
	t195 = t200 * t226 - t220 * t231;
	t194 = t200 * t224 + t220 * t230;
	t193 = t199 * t224 + t217 * t230;
	t1 = [t196 * t213 + t198 * t212, -t196 * t212 + t198 * t213, -t193, t196 * t207 + t193 * t228 + t198 * t234 + (t220 * t205 + t206 * t233) * t227 + (-t205 * t233 + t220 * t206) * t225 + t217 * t232 + t220 * pkin(1) + 0; t195 * t213 + t197 * t212, -t195 * t212 + t197 * t213, t194, t195 * t207 - t194 * t228 + t197 * t234 + (t217 * t205 - t206 * t229) * t227 + (t205 * t229 + t217 * t206) * t225 - t220 * t232 + t217 * pkin(1) + 0; t202 * t213 + t203 * t212, -t202 * t212 + t203 * t213, t201, t203 * t234 - t201 * t228 + t202 * t207 + t221 * t222 + qJ(1) + 0 + (t205 * t225 - t206 * t227) * t218; 0, 0, 0, 1;];
	Tc_mdh = t1;
end