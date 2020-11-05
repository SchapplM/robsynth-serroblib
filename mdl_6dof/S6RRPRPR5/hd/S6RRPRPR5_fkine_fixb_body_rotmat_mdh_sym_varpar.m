% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t120 = cos(qJ(1));
	t119 = sin(qJ(1));
	t1 = [t120, -t119, 0, 0; t119, t120, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t121 = sin(pkin(6));
	t124 = sin(qJ(1));
	t132 = t124 * t121;
	t123 = sin(qJ(2));
	t131 = t124 * t123;
	t125 = cos(qJ(2));
	t130 = t124 * t125;
	t126 = cos(qJ(1));
	t129 = t126 * t121;
	t128 = t126 * t123;
	t127 = t126 * t125;
	t122 = cos(pkin(6));
	t1 = [-t122 * t131 + t127, -t122 * t130 - t128, t132, t126 * pkin(1) + pkin(8) * t132 + 0; t122 * t128 + t130, t122 * t127 - t131, -t129, t124 * pkin(1) - pkin(8) * t129 + 0; t121 * t123, t121 * t125, t122, t122 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t152 = pkin(2) * sin(qJ(2));
	t145 = qJ(2) + pkin(11);
	t151 = cos(qJ(1));
	t150 = sin(qJ(1));
	t148 = pkin(8) + qJ(3);
	t147 = cos(pkin(6));
	t146 = sin(pkin(6));
	t144 = cos(t145);
	t143 = sin(t145);
	t142 = pkin(6) - t145;
	t141 = pkin(6) + t145;
	t140 = cos(qJ(2)) * pkin(2) + pkin(1);
	t139 = cos(t142);
	t138 = cos(t141);
	t137 = sin(t142);
	t136 = sin(t141);
	t135 = t138 + t139;
	t134 = -t136 + t137;
	t133 = -t146 * t148 + t147 * t152;
	t1 = [t151 * t144 + t150 * t134 / 0.2e1, -t151 * t143 - t150 * t135 / 0.2e1, t150 * t146, -t133 * t150 + t151 * t140 + 0; t150 * t144 - t151 * t134 / 0.2e1, -t150 * t143 + t151 * t135 / 0.2e1, -t151 * t146, t151 * t133 + t150 * t140 + 0; t139 / 0.2e1 - t138 / 0.2e1, t137 / 0.2e1 + t136 / 0.2e1, t147, t146 * t152 + t147 * t148 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t166 = sin(pkin(6));
	t170 = sin(qJ(4));
	t181 = t166 * t170;
	t173 = cos(qJ(4));
	t180 = t166 * t173;
	t175 = cos(qJ(1));
	t179 = t166 * t175;
	t164 = qJ(2) + pkin(11);
	t162 = sin(t164);
	t172 = sin(qJ(1));
	t178 = t172 * t162;
	t177 = t175 * t162;
	t165 = sin(pkin(11));
	t167 = cos(pkin(11));
	t158 = t167 * pkin(3) + t165 * pkin(9) + pkin(2);
	t159 = -t165 * pkin(3) + t167 * pkin(9);
	t171 = sin(qJ(2));
	t174 = cos(qJ(2));
	t176 = t158 * t171 - t159 * t174;
	t169 = pkin(8) + qJ(3);
	t168 = cos(pkin(6));
	t163 = cos(t164);
	t161 = pkin(6) - t164;
	t160 = pkin(6) + t164;
	t157 = cos(t160) + cos(t161);
	t156 = t172 * t163 + t168 * t177;
	t155 = -t175 * t163 + t168 * t178;
	t154 = t158 * t174 + t159 * t171 + pkin(1);
	t153 = t166 * t169 - t176 * t168;
	t1 = [-t155 * t173 + t172 * t181, t155 * t170 + t172 * t180, t177 + t172 * t157 / 0.2e1, t153 * t172 + t154 * t175 + 0; t156 * t173 - t170 * t179, -t156 * t170 - t173 * t179, t178 - t175 * t157 / 0.2e1, -t153 * t175 + t154 * t172 + 0; t162 * t180 + t168 * t170, -t162 * t181 + t168 * t173, -sin(t161) / 0.2e1 - sin(t160) / 0.2e1, t176 * t166 + t168 * t169 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->44), mult. (170->73), div. (0->0), fcn. (204->14), ass. (0->36)
	t189 = qJ(2) + pkin(11);
	t188 = cos(t189);
	t192 = sin(pkin(6));
	t216 = t188 * t192;
	t196 = sin(qJ(4));
	t215 = t192 * t196;
	t195 = cos(pkin(6));
	t198 = sin(qJ(1));
	t214 = t195 * t198;
	t201 = cos(qJ(1));
	t213 = t195 * t201;
	t199 = cos(qJ(4));
	t212 = t198 * t199;
	t211 = t199 * t201;
	t190 = sin(pkin(12));
	t210 = t190 * t214;
	t209 = t190 * t215;
	t193 = cos(pkin(12));
	t208 = t193 * t215;
	t207 = t193 * t212;
	t206 = t190 * t211;
	t205 = t193 * t211;
	t191 = sin(pkin(11));
	t194 = cos(pkin(11));
	t203 = pkin(4) * t199 + qJ(5) * t196 + pkin(3);
	t184 = t191 * pkin(9) + t203 * t194 + pkin(2);
	t185 = -t194 * pkin(9) + t203 * t191;
	t197 = sin(qJ(2));
	t200 = cos(qJ(2));
	t204 = t184 * t197 + t185 * t200;
	t202 = t196 * pkin(4) - qJ(5) * t199 + pkin(8) + qJ(3);
	t187 = sin(t189);
	t186 = t192 * t187 * t199 + t195 * t196;
	t183 = t184 * t200 - t185 * t197 + pkin(1);
	t182 = -t192 * t202 + t204 * t195;
	t1 = [(t205 + t210) * t188 + (t190 * t201 - t195 * t207) * t187 + t198 * t208, (t193 * t214 - t206) * t188 + (t193 * t201 + t199 * t210) * t187 - t198 * t209, -(t187 * t214 - t201 * t188) * t196 - t192 * t212, -t182 * t198 + t183 * t201 + 0; (-t190 * t213 + t207) * t188 + (t190 * t198 + t195 * t205) * t187 - t201 * t208, (-t190 * t212 - t193 * t213) * t188 + (t193 * t198 - t195 * t206) * t187 + t201 * t209, (t187 * t213 + t198 * t188) * t196 + t192 * t211, t182 * t201 + t183 * t198 + 0; t186 * t193 - t190 * t216, -t186 * t190 - t193 * t216, t187 * t215 - t195 * t199, t204 * t192 + t202 * t195 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:09:25
	% EndTime: 2020-11-04 22:09:25
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (159->50), mult. (187->69), div. (0->0), fcn. (216->18), ass. (0->43)
	t261 = sin(pkin(12)) * pkin(5);
	t226 = pkin(9) + t261;
	t236 = sin(pkin(11));
	t238 = cos(pkin(11));
	t227 = cos(pkin(12)) * pkin(5) + pkin(4);
	t240 = qJ(5) + pkin(10);
	t241 = sin(qJ(4));
	t244 = cos(qJ(4));
	t249 = t227 * t244 + t241 * t240;
	t262 = t236 * (pkin(3) + t249) - t226 * t238;
	t260 = t238 * pkin(3) + pkin(2);
	t234 = qJ(2) + pkin(11);
	t231 = cos(t234);
	t237 = sin(pkin(6));
	t258 = t231 * t237;
	t257 = t237 * t241;
	t256 = t237 * t244;
	t239 = cos(pkin(6));
	t255 = t239 * t244;
	t229 = sin(t234);
	t243 = sin(qJ(1));
	t254 = t243 * t229;
	t253 = t243 * t231;
	t252 = t243 * t244;
	t246 = cos(qJ(1));
	t251 = t246 * t229;
	t250 = t246 * t231;
	t247 = t227 * t241 - t240 * t244 + pkin(8) + qJ(3);
	t245 = cos(qJ(2));
	t242 = sin(qJ(2));
	t233 = pkin(12) + qJ(6);
	t230 = cos(t233);
	t228 = sin(t233);
	t225 = -t229 * t255 + t257;
	t224 = t229 * t256 + t239 * t241;
	t223 = t239 * t253 + t251;
	t222 = t239 * t250 - t254;
	t221 = t243 * t225 + t244 * t250;
	t220 = -t246 * t225 + t231 * t252;
	t219 = t226 * t236 + t249 * t238 + t260;
	t218 = t219 * t245 - t262 * t242 + pkin(1);
	t217 = -t247 * t237 + (t219 * t242 + t262 * t245) * t239;
	t1 = [t221 * t230 + t223 * t228, -t221 * t228 + t223 * t230, -(t239 * t254 - t250) * t241 - t237 * t252, -t217 * t243 + t218 * t246 + 0; t220 * t230 - t228 * t222, -t220 * t228 - t230 * t222, (t239 * t251 + t253) * t241 + t246 * t256, t217 * t246 + t218 * t243 + 0; t224 * t230 - t228 * t258, -t224 * t228 - t230 * t258, t229 * t257 - t255, pkin(7) + 0 + (-sin(pkin(6) - t234) / 0.2e1 - sin(pkin(6) + t234) / 0.2e1) * t261 + t247 * t239 + (t249 * t229 + (t236 * pkin(9) + t260) * t242 - (-t236 * pkin(3) + t238 * pkin(9)) * t245) * t237; 0, 0, 0, 1;];
	Tc_mdh = t1;
end