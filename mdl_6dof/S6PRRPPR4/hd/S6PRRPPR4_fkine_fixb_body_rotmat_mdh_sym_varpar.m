% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t121 = cos(pkin(10));
	t120 = sin(pkin(10));
	t1 = [t121, -t120, 0, 0; t120, t121, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t122 = sin(pkin(10));
	t123 = sin(pkin(6));
	t131 = t122 * t123;
	t124 = cos(pkin(10));
	t130 = t124 * t123;
	t125 = cos(pkin(6));
	t126 = sin(qJ(2));
	t129 = t125 * t126;
	t127 = cos(qJ(2));
	t128 = t125 * t127;
	t1 = [-t122 * t129 + t124 * t127, -t122 * t128 - t124 * t126, t131, t124 * pkin(1) + pkin(7) * t131 + 0; t122 * t127 + t124 * t129, -t122 * t126 + t124 * t128, -t130, t122 * pkin(1) - pkin(7) * t130 + 0; t123 * t126, t123 * t127, t125, t125 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t135 = sin(pkin(6));
	t148 = t135 * pkin(7);
	t134 = sin(pkin(10));
	t137 = cos(pkin(6));
	t147 = t134 * t137;
	t138 = sin(qJ(3));
	t146 = t135 * t138;
	t140 = cos(qJ(3));
	t145 = t135 * t140;
	t136 = cos(pkin(10));
	t144 = t136 * t137;
	t139 = sin(qJ(2));
	t143 = t137 * t139;
	t141 = cos(qJ(2));
	t142 = t137 * t141;
	t133 = -t134 * t143 + t136 * t141;
	t132 = t134 * t141 + t136 * t143;
	t1 = [t133 * t140 + t134 * t146, -t133 * t138 + t134 * t145, t134 * t142 + t136 * t139, (t136 * pkin(2) + pkin(8) * t147) * t141 + (-pkin(2) * t147 + t136 * pkin(8)) * t139 + t134 * t148 + t136 * pkin(1) + 0; t132 * t140 - t136 * t146, -t132 * t138 - t136 * t145, t134 * t139 - t136 * t142, (t134 * pkin(2) - pkin(8) * t144) * t141 + (pkin(2) * t144 + t134 * pkin(8)) * t139 - t136 * t148 + t134 * pkin(1) + 0; t137 * t138 + t139 * t145, t137 * t140 - t139 * t146, -t135 * t141, t137 * pkin(7) + qJ(1) + 0 + (pkin(2) * t139 - pkin(8) * t141) * t135; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->39), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->29)
	t161 = sin(pkin(6));
	t176 = t161 * pkin(7);
	t160 = sin(pkin(10));
	t164 = cos(pkin(6));
	t175 = t160 * t164;
	t165 = sin(qJ(3));
	t174 = t161 * t165;
	t167 = cos(qJ(3));
	t173 = t161 * t167;
	t168 = cos(qJ(2));
	t172 = t161 * t168;
	t163 = cos(pkin(10));
	t171 = t163 * t164;
	t166 = sin(qJ(2));
	t170 = t164 * t166;
	t169 = t164 * t168;
	t162 = cos(pkin(11));
	t159 = sin(pkin(11));
	t158 = t164 * t165 + t166 * t173;
	t157 = -t164 * t167 + t166 * t174;
	t156 = -t160 * t170 + t163 * t168;
	t155 = t160 * t169 + t163 * t166;
	t154 = t160 * t168 + t163 * t170;
	t153 = t160 * t166 - t163 * t169;
	t152 = -t156 * t165 + t160 * t173;
	t151 = t156 * t167 + t160 * t174;
	t150 = t154 * t167 - t163 * t174;
	t149 = t154 * t165 + t163 * t173;
	t1 = [t151 * t162 + t155 * t159, -t151 * t159 + t155 * t162, -t152, t151 * pkin(3) - t152 * qJ(4) + (t163 * pkin(2) + pkin(8) * t175) * t168 + (-pkin(2) * t175 + t163 * pkin(8)) * t166 + t160 * t176 + t163 * pkin(1) + 0; t150 * t162 + t153 * t159, -t150 * t159 + t153 * t162, t149, t150 * pkin(3) + t149 * qJ(4) + (t160 * pkin(2) - pkin(8) * t171) * t168 + (pkin(2) * t171 + t160 * pkin(8)) * t166 - t163 * t176 + t160 * pkin(1) + 0; t158 * t162 - t159 * t172, -t158 * t159 - t162 * t172, t157, t158 * pkin(3) + t164 * pkin(7) + t157 * qJ(4) + qJ(1) + 0 + (pkin(2) * t166 - pkin(8) * t168) * t161; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:37
	% EndTime: 2020-11-04 21:07:38
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (83->41), mult. (170->70), div. (0->0), fcn. (204->10), ass. (0->32)
	t184 = sin(pkin(6));
	t187 = cos(pkin(6));
	t182 = sin(pkin(11));
	t185 = cos(pkin(11));
	t181 = -t182 * pkin(4) + qJ(5) * t185 - pkin(8);
	t189 = sin(qJ(2));
	t191 = cos(qJ(2));
	t180 = t185 * pkin(4) + qJ(5) * t182 + pkin(3);
	t188 = sin(qJ(3));
	t190 = cos(qJ(3));
	t195 = qJ(4) * t188 + t180 * t190 + pkin(2);
	t193 = t181 * t191 + t195 * t189;
	t194 = qJ(4) * t190 - t180 * t188 - pkin(7);
	t210 = t194 * t184 + t193 * t187;
	t206 = t184 * t188;
	t205 = t184 * t190;
	t204 = t184 * t191;
	t203 = t187 * t189;
	t202 = t187 * t190;
	t201 = t187 * t191;
	t200 = t190 * t191;
	t183 = sin(pkin(10));
	t199 = t183 * t202;
	t186 = cos(pkin(10));
	t198 = t186 * t202;
	t197 = t183 * t200;
	t196 = t186 * t200;
	t192 = -t181 * t189 + t195 * t191 + pkin(1);
	t179 = t187 * t188 + t189 * t205;
	t178 = t182 * t201 + t185 * t206;
	t177 = t182 * t206 - t185 * t201;
	t1 = [(t182 * t186 - t185 * t199) * t189 + t185 * t196 + t183 * t178, -(t183 * t203 - t186 * t191) * t188 - t183 * t205, (-t182 * t199 - t186 * t185) * t189 + t182 * t196 + t183 * t177, -t210 * t183 + t192 * t186 + 0; (t183 * t182 + t185 * t198) * t189 + t185 * t197 - t186 * t178, (t183 * t191 + t186 * t203) * t188 + t186 * t205, (t182 * t198 - t183 * t185) * t189 + t182 * t197 - t186 * t177, t192 * t183 + t210 * t186 + 0; t179 * t185 - t182 * t204, t189 * t206 - t202, t179 * t182 + t185 * t204, t193 * t184 - t194 * t187 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:07:38
	% EndTime: 2020-11-04 21:07:38
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (130->50), mult. (254->87), div. (0->0), fcn. (314->12), ass. (0->41)
	t222 = sin(pkin(6));
	t225 = cos(pkin(6));
	t220 = sin(pkin(11));
	t223 = cos(pkin(11));
	t233 = pkin(4) + pkin(5);
	t219 = qJ(5) * t223 - t233 * t220 - pkin(8);
	t229 = sin(qJ(2));
	t232 = cos(qJ(2));
	t218 = qJ(5) * t220 + t233 * t223 + pkin(3);
	t226 = qJ(4) - pkin(9);
	t228 = sin(qJ(3));
	t231 = cos(qJ(3));
	t243 = t218 * t231 + t226 * t228 + pkin(2);
	t235 = t219 * t232 + t243 * t229;
	t242 = t218 * t228 - t226 * t231 + pkin(7);
	t258 = t242 * t222 - t235 * t225;
	t221 = sin(pkin(10));
	t244 = t229 * t231;
	t249 = t222 * t228;
	t238 = t225 * t244 - t249;
	t224 = cos(pkin(10));
	t245 = t229 * t224;
	t246 = t224 * t231;
	t255 = t238 * t221 * t220 + (t223 * t221 * t225 - t220 * t246) * t232 + t223 * t245;
	t252 = t220 * t232;
	t251 = t221 * t229;
	t250 = t221 * t231;
	t248 = t223 * t232;
	t247 = t224 * t225;
	t241 = (-t220 * t247 + t223 * t250) * t232 + t220 * t251;
	t239 = -(t220 * t250 + t223 * t247) * t232 + t223 * t251;
	t227 = sin(qJ(6));
	t237 = t238 * t227;
	t230 = cos(qJ(6));
	t236 = t238 * t230;
	t234 = -t219 * t229 + t243 * t232 + pkin(1);
	t217 = t222 * t244 + t225 * t228;
	t213 = t217 * t223 - t222 * t252;
	t212 = t217 * t220 + t222 * t248;
	t211 = (t220 * t229 + t231 * t248) * t224 + (-t238 * t223 + t225 * t252) * t221;
	t1 = [t211 * t230 - t255 * t227, -t211 * t227 - t255 * t230, (-t224 * t232 + t225 * t251) * t228 + t222 * t250, t258 * t221 + t234 * t224 + 0; t241 * t230 - t239 * t227 + (t220 * t237 + t223 * t236) * t224, -t239 * t230 - t241 * t227 + (t220 * t236 - t223 * t237) * t224, -(t221 * t232 + t225 * t245) * t228 - t222 * t246, t234 * t221 - t258 * t224 + 0; t212 * t227 + t213 * t230, t212 * t230 - t213 * t227, t225 * t231 - t229 * t249, t235 * t222 + t242 * t225 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end