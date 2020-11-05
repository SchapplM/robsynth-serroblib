% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t115 = cos(qJ(1));
	t114 = sin(qJ(1));
	t1 = [t115, -t114, 0, 0; t114, t115, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t116 = sin(pkin(6));
	t119 = sin(qJ(1));
	t127 = t119 * t116;
	t118 = sin(qJ(2));
	t126 = t119 * t118;
	t120 = cos(qJ(2));
	t125 = t119 * t120;
	t121 = cos(qJ(1));
	t124 = t121 * t116;
	t123 = t121 * t118;
	t122 = t121 * t120;
	t117 = cos(pkin(6));
	t1 = [-t117 * t126 + t122, -t117 * t125 - t123, t127, t121 * pkin(1) + pkin(8) * t127 + 0; t117 * t123 + t125, t117 * t122 - t126, -t124, t119 * pkin(1) - pkin(8) * t124 + 0; t116 * t118, t116 * t120, t117, t117 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t147 = pkin(2) * sin(qJ(2));
	t140 = qJ(2) + pkin(11);
	t146 = cos(qJ(1));
	t145 = sin(qJ(1));
	t143 = pkin(8) + qJ(3);
	t142 = cos(pkin(6));
	t141 = sin(pkin(6));
	t139 = cos(t140);
	t138 = sin(t140);
	t137 = pkin(6) - t140;
	t136 = pkin(6) + t140;
	t135 = cos(qJ(2)) * pkin(2) + pkin(1);
	t134 = cos(t137);
	t133 = cos(t136);
	t132 = sin(t137);
	t131 = sin(t136);
	t130 = t133 + t134;
	t129 = -t131 + t132;
	t128 = -t141 * t143 + t142 * t147;
	t1 = [t146 * t139 + t145 * t129 / 0.2e1, -t146 * t138 - t145 * t130 / 0.2e1, t145 * t141, -t128 * t145 + t146 * t135 + 0; t145 * t139 - t146 * t129 / 0.2e1, -t145 * t138 + t146 * t130 / 0.2e1, -t146 * t141, t146 * t128 + t145 * t135 + 0; t134 / 0.2e1 - t133 / 0.2e1, t132 / 0.2e1 + t131 / 0.2e1, t142, t141 * t147 + t142 * t143 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t161 = sin(pkin(6));
	t165 = sin(qJ(4));
	t176 = t161 * t165;
	t168 = cos(qJ(4));
	t175 = t161 * t168;
	t170 = cos(qJ(1));
	t174 = t161 * t170;
	t159 = qJ(2) + pkin(11);
	t157 = sin(t159);
	t167 = sin(qJ(1));
	t173 = t167 * t157;
	t172 = t170 * t157;
	t160 = sin(pkin(11));
	t162 = cos(pkin(11));
	t153 = t162 * pkin(3) + t160 * pkin(9) + pkin(2);
	t154 = -t160 * pkin(3) + t162 * pkin(9);
	t166 = sin(qJ(2));
	t169 = cos(qJ(2));
	t171 = t153 * t166 - t154 * t169;
	t164 = pkin(8) + qJ(3);
	t163 = cos(pkin(6));
	t158 = cos(t159);
	t156 = pkin(6) - t159;
	t155 = pkin(6) + t159;
	t152 = cos(t155) + cos(t156);
	t151 = t167 * t158 + t163 * t172;
	t150 = -t170 * t158 + t163 * t173;
	t149 = t153 * t169 + t154 * t166 + pkin(1);
	t148 = t161 * t164 - t171 * t163;
	t1 = [-t150 * t168 + t167 * t176, t150 * t165 + t167 * t175, t172 + t167 * t152 / 0.2e1, t148 * t167 + t149 * t170 + 0; t151 * t168 - t165 * t174, -t151 * t165 - t168 * t174, t173 - t170 * t152 / 0.2e1, -t148 * t170 + t149 * t167 + 0; t157 * t175 + t163 * t165, -t157 * t176 + t163 * t168, -sin(t156) / 0.2e1 - sin(t155) / 0.2e1, t171 * t161 + t163 * t164 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:10
	% EndTime: 2020-11-04 22:14:11
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (102->44), mult. (170->73), div. (0->0), fcn. (204->14), ass. (0->39)
	t184 = qJ(2) + pkin(11);
	t182 = sin(t184);
	t188 = cos(pkin(6));
	t214 = t182 * t188;
	t186 = sin(pkin(6));
	t190 = sin(qJ(4));
	t213 = t186 * t190;
	t193 = cos(qJ(5));
	t212 = t186 * t193;
	t189 = sin(qJ(5));
	t192 = sin(qJ(1));
	t211 = t189 * t192;
	t196 = cos(qJ(1));
	t210 = t189 * t196;
	t194 = cos(qJ(4));
	t209 = t192 * t194;
	t208 = t193 * t192;
	t207 = t193 * t196;
	t206 = t194 * t196;
	t205 = t189 * t213;
	t204 = t190 * t212;
	t203 = t194 * t208;
	t202 = t189 * t209;
	t201 = t189 * t206;
	t200 = t193 * t206;
	t185 = sin(pkin(11));
	t187 = cos(pkin(11));
	t198 = pkin(4) * t194 + pkin(10) * t190 + pkin(3);
	t179 = t185 * pkin(9) + t198 * t187 + pkin(2);
	t180 = t187 * pkin(9) - t198 * t185;
	t191 = sin(qJ(2));
	t195 = cos(qJ(2));
	t199 = t179 * t191 - t180 * t195;
	t197 = t190 * pkin(4) - pkin(10) * t194 + pkin(8) + qJ(3);
	t183 = cos(t184);
	t181 = t186 * t182 * t194 + t188 * t190;
	t178 = t179 * t195 + t180 * t191 + pkin(1);
	t177 = t186 * t197 - t199 * t188;
	t1 = [(t188 * t211 + t200) * t183 + (-t188 * t203 + t210) * t182 + t192 * t204, (t188 * t208 - t201) * t183 + (t188 * t202 + t207) * t182 - t192 * t205, -(-t196 * t183 + t192 * t214) * t190 - t186 * t209, t177 * t192 + t178 * t196 + 0; (-t188 * t210 + t203) * t183 + (t188 * t200 + t211) * t182 - t196 * t204, (-t188 * t207 - t202) * t183 + (-t188 * t201 + t208) * t182 + t196 * t205, (t192 * t183 + t196 * t214) * t190 + t186 * t206, -t177 * t196 + t178 * t192 + 0; -t186 * t183 * t189 + t181 * t193, -t181 * t189 - t183 * t212, t182 * t213 - t188 * t194, t199 * t186 + t188 * t197 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:11
	% EndTime: 2020-11-04 22:14:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (162->48), mult. (169->65), div. (0->0), fcn. (212->18), ass. (0->46)
	t236 = qJ(2) + pkin(11);
	t232 = pkin(6) + t236;
	t233 = pkin(6) - t236;
	t228 = cos(t232) + cos(t233);
	t250 = cos(qJ(1));
	t234 = sin(t236);
	t246 = sin(qJ(1));
	t253 = t246 * t234;
	t221 = t253 - t250 * t228 / 0.2e1;
	t243 = sin(qJ(5));
	t259 = t221 * t243;
	t252 = t250 * t234;
	t222 = t252 + t246 * t228 / 0.2e1;
	t258 = t222 * t243;
	t227 = -sin(t233) / 0.2e1 - sin(t232) / 0.2e1;
	t257 = t227 * t243;
	t238 = sin(pkin(6));
	t244 = sin(qJ(4));
	t256 = t238 * t244;
	t248 = cos(qJ(4));
	t255 = t238 * t248;
	t254 = t238 * t250;
	t237 = sin(pkin(11));
	t239 = cos(pkin(11));
	t229 = t239 * pkin(3) + t237 * pkin(9) + pkin(2);
	t230 = -t237 * pkin(3) + t239 * pkin(9);
	t245 = sin(qJ(2));
	t249 = cos(qJ(2));
	t251 = t229 * t245 - t230 * t249;
	t247 = cos(qJ(5));
	t242 = pkin(8) + qJ(3);
	t241 = -qJ(6) - pkin(10);
	t240 = cos(pkin(6));
	t235 = cos(t236);
	t231 = t247 * pkin(5) + pkin(4);
	t226 = t234 * t255 + t240 * t244;
	t225 = t234 * t256 - t240 * t248;
	t224 = t246 * t235 + t240 * t252;
	t223 = -t250 * t235 + t240 * t253;
	t220 = t229 * t249 + t230 * t245 + pkin(1);
	t219 = -t223 * t248 + t246 * t256;
	t218 = t224 * t248 - t244 * t254;
	t217 = t224 * t244 + t248 * t254;
	t216 = t223 * t244 + t246 * t255;
	t215 = t238 * t242 - t251 * t240;
	t1 = [t219 * t247 + t258, -t219 * t243 + t222 * t247, -t216, pkin(5) * t258 + t215 * t246 + t216 * t241 + t219 * t231 + t220 * t250 + 0; t218 * t247 + t259, -t218 * t243 + t221 * t247, t217, pkin(5) * t259 - t215 * t250 - t217 * t241 + t218 * t231 + t220 * t246 + 0; t226 * t247 + t257, -t226 * t243 + t227 * t247, t225, pkin(5) * t257 - t225 * t241 + t226 * t231 + t251 * t238 + t240 * t242 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end