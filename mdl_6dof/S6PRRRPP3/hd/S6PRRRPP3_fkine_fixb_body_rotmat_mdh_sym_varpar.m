% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t123 = cos(pkin(10));
	t122 = sin(pkin(10));
	t1 = [t123, -t122, 0, 0; t122, t123, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t124 = sin(pkin(10));
	t125 = sin(pkin(6));
	t133 = t124 * t125;
	t126 = cos(pkin(10));
	t132 = t126 * t125;
	t127 = cos(pkin(6));
	t128 = sin(qJ(2));
	t131 = t127 * t128;
	t129 = cos(qJ(2));
	t130 = t127 * t129;
	t1 = [-t124 * t131 + t126 * t129, -t124 * t130 - t126 * t128, t133, t126 * pkin(1) + pkin(7) * t133 + 0; t124 * t129 + t126 * t131, -t124 * t128 + t126 * t130, -t132, t124 * pkin(1) - pkin(7) * t132 + 0; t125 * t128, t125 * t129, t127, t127 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t137 = sin(pkin(6));
	t150 = t137 * pkin(7);
	t136 = sin(pkin(10));
	t139 = cos(pkin(6));
	t149 = t136 * t139;
	t140 = sin(qJ(3));
	t148 = t137 * t140;
	t142 = cos(qJ(3));
	t147 = t137 * t142;
	t138 = cos(pkin(10));
	t146 = t138 * t139;
	t141 = sin(qJ(2));
	t145 = t139 * t141;
	t143 = cos(qJ(2));
	t144 = t139 * t143;
	t135 = -t136 * t145 + t138 * t143;
	t134 = t136 * t143 + t138 * t145;
	t1 = [t135 * t142 + t136 * t148, -t135 * t140 + t136 * t147, t136 * t144 + t138 * t141, (t138 * pkin(2) + pkin(8) * t149) * t143 + (-pkin(2) * t149 + t138 * pkin(8)) * t141 + t136 * t150 + t138 * pkin(1) + 0; t134 * t142 - t138 * t148, -t134 * t140 - t138 * t147, t136 * t141 - t138 * t144, (t136 * pkin(2) - pkin(8) * t146) * t143 + (pkin(2) * t146 + t136 * pkin(8)) * t141 - t138 * t150 + t136 * pkin(1) + 0; t139 * t140 + t141 * t147, t139 * t142 - t141 * t148, -t137 * t143, t139 * pkin(7) + qJ(1) + 0 + (pkin(2) * t141 - pkin(8) * t143) * t137; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t158 = sin(pkin(6));
	t160 = cos(pkin(6));
	t163 = sin(qJ(2));
	t166 = cos(qJ(2));
	t162 = sin(qJ(3));
	t165 = cos(qJ(3));
	t170 = t165 * pkin(3) + t162 * pkin(9) + pkin(2);
	t168 = -pkin(8) * t166 + t170 * t163;
	t169 = t162 * pkin(3) - t165 * pkin(9) + pkin(7);
	t181 = t169 * t158 - t168 * t160;
	t177 = t158 * t165;
	t176 = t158 * t166;
	t175 = t160 * t163;
	t174 = t160 * t166;
	t173 = t162 * t158;
	t172 = t163 * t165;
	t171 = t165 * t166;
	t167 = pkin(8) * t163 + t170 * t166 + pkin(1);
	t164 = cos(qJ(4));
	t161 = sin(qJ(4));
	t159 = cos(pkin(10));
	t157 = sin(pkin(10));
	t156 = t160 * t172 - t173;
	t155 = t158 * t172 + t160 * t162;
	t154 = t157 * t174 + t159 * t163;
	t153 = t157 * t163 - t159 * t174;
	t152 = -t157 * t156 + t159 * t171;
	t151 = t159 * t156 + t157 * t171;
	t1 = [t152 * t164 + t154 * t161, -t152 * t161 + t154 * t164, -(t157 * t175 - t159 * t166) * t162 - t157 * t177, t181 * t157 + t167 * t159 + 0; t151 * t164 + t153 * t161, -t151 * t161 + t153 * t164, (t157 * t166 + t159 * t175) * t162 + t159 * t177, t167 * t157 - t181 * t159 + 0; t155 * t164 - t161 * t176, -t155 * t161 - t164 * t176, -t160 * t165 + t163 * t173, t168 * t158 + t169 * t160 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:15
	% EndTime: 2020-11-04 21:14:16
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (83->36), mult. (180->60), div. (0->0), fcn. (214->10), ass. (0->31)
	t192 = sin(pkin(6));
	t194 = cos(pkin(6));
	t197 = sin(qJ(2));
	t200 = cos(qJ(2));
	t195 = sin(qJ(4));
	t198 = cos(qJ(4));
	t203 = t195 * pkin(4) - qJ(5) * t198 + pkin(8);
	t196 = sin(qJ(3));
	t199 = cos(qJ(3));
	t216 = pkin(4) * t198 + qJ(5) * t195 + pkin(3);
	t204 = t196 * pkin(9) + t199 * t216 + pkin(2);
	t202 = t204 * t197 - t203 * t200;
	t215 = -t199 * pkin(9) + t196 * t216 + pkin(7);
	t218 = t215 * t192 - t202 * t194;
	t213 = t192 * t199;
	t212 = t192 * t200;
	t211 = t194 * t197;
	t210 = t194 * t200;
	t209 = t196 * t192;
	t208 = t197 * t199;
	t207 = t199 * t200;
	t201 = t203 * t197 + t204 * t200 + pkin(1);
	t193 = cos(pkin(10));
	t191 = sin(pkin(10));
	t187 = t194 * t208 - t209;
	t186 = t192 * t208 + t194 * t196;
	t185 = t191 * t210 + t193 * t197;
	t184 = t191 * t197 - t193 * t210;
	t183 = -t191 * t187 + t193 * t207;
	t182 = t193 * t187 + t191 * t207;
	t1 = [-(t191 * t211 - t193 * t200) * t196 - t191 * t213, -t183 * t198 - t185 * t195, t183 * t195 - t185 * t198, t218 * t191 + t201 * t193 + 0; (t191 * t200 + t193 * t211) * t196 + t193 * t213, -t182 * t198 - t184 * t195, t182 * t195 - t184 * t198, t201 * t191 - t218 * t193 + 0; -t194 * t199 + t197 * t209, -t186 * t198 + t195 * t212, t186 * t195 + t198 * t212, t202 * t192 + t215 * t194 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:14:16
	% EndTime: 2020-11-04 21:14:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (104->38), mult. (182->60), div. (0->0), fcn. (216->10), ass. (0->33)
	t229 = sin(pkin(6));
	t231 = cos(pkin(6));
	t235 = sin(qJ(2));
	t238 = cos(qJ(2));
	t232 = qJ(6) + pkin(4);
	t233 = sin(qJ(4));
	t236 = cos(qJ(4));
	t242 = qJ(5) * t236 - t232 * t233 - pkin(8);
	t234 = sin(qJ(3));
	t237 = cos(qJ(3));
	t239 = pkin(5) + pkin(9);
	t255 = qJ(5) * t233 + t232 * t236 + pkin(3);
	t243 = t239 * t234 + t237 * t255 + pkin(2);
	t241 = t243 * t235 + t242 * t238;
	t254 = t234 * t255 - t239 * t237 + pkin(7);
	t257 = t254 * t229 - t241 * t231;
	t253 = t229 * t237;
	t252 = t229 * t238;
	t251 = t231 * t235;
	t250 = t231 * t238;
	t248 = t234 * t229;
	t247 = t235 * t237;
	t246 = t237 * t238;
	t240 = -t242 * t235 + t243 * t238 + pkin(1);
	t230 = cos(pkin(10));
	t228 = sin(pkin(10));
	t224 = t231 * t247 - t248;
	t223 = t229 * t247 + t231 * t234;
	t222 = t228 * t250 + t230 * t235;
	t221 = t228 * t235 - t230 * t250;
	t220 = -t228 * t224 + t230 * t246;
	t219 = t230 * t224 + t228 * t246;
	t1 = [-(t228 * t251 - t230 * t238) * t234 - t228 * t253, t220 * t233 - t222 * t236, t220 * t236 + t222 * t233, t257 * t228 + t240 * t230 + 0; (t228 * t238 + t230 * t251) * t234 + t230 * t253, t219 * t233 - t221 * t236, t219 * t236 + t221 * t233, t240 * t228 - t257 * t230 + 0; -t231 * t237 + t235 * t248, t223 * t233 + t236 * t252, t223 * t236 - t233 * t252, t241 * t229 + t254 * t231 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end