% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPP2 (for one body)
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
% Datum: 2020-11-04 21:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:52
	% EndTime: 2020-11-04 21:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:52
	% EndTime: 2020-11-04 21:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t124 = cos(pkin(10));
	t123 = sin(pkin(10));
	t1 = [t124, -t123, 0, 0; t123, t124, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:52
	% EndTime: 2020-11-04 21:13:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t125 = sin(pkin(10));
	t126 = sin(pkin(6));
	t134 = t125 * t126;
	t127 = cos(pkin(10));
	t133 = t127 * t126;
	t128 = cos(pkin(6));
	t129 = sin(qJ(2));
	t132 = t128 * t129;
	t130 = cos(qJ(2));
	t131 = t128 * t130;
	t1 = [-t125 * t132 + t127 * t130, -t125 * t131 - t127 * t129, t134, t127 * pkin(1) + pkin(7) * t134 + 0; t125 * t130 + t127 * t132, -t125 * t129 + t127 * t131, -t133, t125 * pkin(1) - pkin(7) * t133 + 0; t126 * t129, t126 * t130, t128, t128 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:52
	% EndTime: 2020-11-04 21:13:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t138 = sin(pkin(6));
	t151 = t138 * pkin(7);
	t137 = sin(pkin(10));
	t140 = cos(pkin(6));
	t150 = t137 * t140;
	t141 = sin(qJ(3));
	t149 = t138 * t141;
	t143 = cos(qJ(3));
	t148 = t138 * t143;
	t139 = cos(pkin(10));
	t147 = t139 * t140;
	t142 = sin(qJ(2));
	t146 = t140 * t142;
	t144 = cos(qJ(2));
	t145 = t140 * t144;
	t136 = -t137 * t146 + t139 * t144;
	t135 = t137 * t144 + t139 * t146;
	t1 = [t136 * t143 + t137 * t149, -t136 * t141 + t137 * t148, t137 * t145 + t139 * t142, (t139 * pkin(2) + pkin(8) * t150) * t144 + (-pkin(2) * t150 + t139 * pkin(8)) * t142 + t137 * t151 + t139 * pkin(1) + 0; t135 * t143 - t139 * t149, -t135 * t141 - t139 * t148, t137 * t142 - t139 * t145, (t137 * pkin(2) - pkin(8) * t147) * t144 + (pkin(2) * t147 + t137 * pkin(8)) * t142 - t139 * t151 + t137 * pkin(1) + 0; t140 * t141 + t142 * t148, t140 * t143 - t142 * t149, -t138 * t144, t140 * pkin(7) + qJ(1) + 0 + (pkin(2) * t142 - pkin(8) * t144) * t138; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:52
	% EndTime: 2020-11-04 21:13:53
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t159 = sin(pkin(6));
	t161 = cos(pkin(6));
	t164 = sin(qJ(2));
	t167 = cos(qJ(2));
	t163 = sin(qJ(3));
	t166 = cos(qJ(3));
	t171 = t166 * pkin(3) + t163 * pkin(9) + pkin(2);
	t169 = -pkin(8) * t167 + t171 * t164;
	t170 = t163 * pkin(3) - t166 * pkin(9) + pkin(7);
	t182 = t170 * t159 - t169 * t161;
	t178 = t159 * t166;
	t177 = t159 * t167;
	t176 = t161 * t164;
	t175 = t161 * t167;
	t174 = t163 * t159;
	t173 = t164 * t166;
	t172 = t166 * t167;
	t168 = pkin(8) * t164 + t171 * t167 + pkin(1);
	t165 = cos(qJ(4));
	t162 = sin(qJ(4));
	t160 = cos(pkin(10));
	t158 = sin(pkin(10));
	t157 = t161 * t173 - t174;
	t156 = t159 * t173 + t161 * t163;
	t155 = t158 * t175 + t160 * t164;
	t154 = t158 * t164 - t160 * t175;
	t153 = -t158 * t157 + t160 * t172;
	t152 = t160 * t157 + t158 * t172;
	t1 = [t153 * t165 + t155 * t162, -t153 * t162 + t155 * t165, -(t158 * t176 - t160 * t167) * t163 - t158 * t178, t182 * t158 + t168 * t160 + 0; t152 * t165 + t154 * t162, -t152 * t162 + t154 * t165, (t158 * t167 + t160 * t176) * t163 + t160 * t178, t168 * t158 - t182 * t160 + 0; t156 * t165 - t162 * t177, -t156 * t162 - t165 * t177, -t161 * t166 + t164 * t174, t169 * t159 + t170 * t161 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:53
	% EndTime: 2020-11-04 21:13:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (83->36), mult. (180->60), div. (0->0), fcn. (214->10), ass. (0->31)
	t193 = sin(pkin(6));
	t195 = cos(pkin(6));
	t198 = sin(qJ(2));
	t201 = cos(qJ(2));
	t196 = sin(qJ(4));
	t199 = cos(qJ(4));
	t204 = t196 * pkin(4) - qJ(5) * t199 + pkin(8);
	t197 = sin(qJ(3));
	t200 = cos(qJ(3));
	t217 = pkin(4) * t199 + qJ(5) * t196 + pkin(3);
	t205 = t197 * pkin(9) + t200 * t217 + pkin(2);
	t203 = t205 * t198 - t204 * t201;
	t216 = -t200 * pkin(9) + t197 * t217 + pkin(7);
	t219 = t216 * t193 - t203 * t195;
	t214 = t193 * t200;
	t213 = t193 * t201;
	t212 = t195 * t198;
	t211 = t195 * t201;
	t210 = t197 * t193;
	t209 = t198 * t200;
	t208 = t200 * t201;
	t202 = t204 * t198 + t205 * t201 + pkin(1);
	t194 = cos(pkin(10));
	t192 = sin(pkin(10));
	t188 = t195 * t209 - t210;
	t187 = t193 * t209 + t195 * t197;
	t186 = t192 * t211 + t194 * t198;
	t185 = t192 * t198 - t194 * t211;
	t184 = -t192 * t188 + t194 * t208;
	t183 = t194 * t188 + t192 * t208;
	t1 = [t184 * t199 + t186 * t196, -(t192 * t212 - t194 * t201) * t197 - t192 * t214, t184 * t196 - t186 * t199, t219 * t192 + t202 * t194 + 0; t183 * t199 + t185 * t196, (t192 * t201 + t194 * t212) * t197 + t194 * t214, t183 * t196 - t185 * t199, t202 * t192 - t219 * t194 + 0; t187 * t199 - t196 * t213, -t195 * t200 + t198 * t210, t187 * t196 + t199 * t213, t203 * t193 + t216 * t195 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:53
	% EndTime: 2020-11-04 21:13:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (104->38), mult. (182->60), div. (0->0), fcn. (216->10), ass. (0->33)
	t230 = sin(pkin(6));
	t232 = cos(pkin(6));
	t236 = sin(qJ(2));
	t239 = cos(qJ(2));
	t234 = sin(qJ(4));
	t237 = cos(qJ(4));
	t240 = pkin(4) + pkin(5);
	t243 = qJ(5) * t237 - t240 * t234 - pkin(8);
	t233 = qJ(6) - pkin(9);
	t235 = sin(qJ(3));
	t238 = cos(qJ(3));
	t256 = qJ(5) * t234 + t240 * t237 + pkin(3);
	t244 = -t233 * t235 + t238 * t256 + pkin(2);
	t242 = t244 * t236 + t243 * t239;
	t255 = t233 * t238 + t235 * t256 + pkin(7);
	t258 = t255 * t230 - t242 * t232;
	t254 = t230 * t238;
	t253 = t230 * t239;
	t252 = t232 * t236;
	t251 = t232 * t239;
	t250 = t235 * t230;
	t249 = t236 * t238;
	t248 = t238 * t239;
	t241 = -t243 * t236 + t244 * t239 + pkin(1);
	t231 = cos(pkin(10));
	t229 = sin(pkin(10));
	t225 = t232 * t249 - t250;
	t224 = t230 * t249 + t232 * t235;
	t223 = t229 * t251 + t231 * t236;
	t222 = t229 * t236 - t231 * t251;
	t221 = -t229 * t225 + t231 * t248;
	t220 = t231 * t225 + t229 * t248;
	t1 = [t221 * t237 + t223 * t234, t221 * t234 - t223 * t237, (t229 * t252 - t231 * t239) * t235 + t229 * t254, t258 * t229 + t241 * t231 + 0; t220 * t237 + t222 * t234, t220 * t234 - t222 * t237, -(t229 * t239 + t231 * t252) * t235 - t231 * t254, t241 * t229 - t258 * t231 + 0; t224 * t237 - t234 * t253, t224 * t234 + t237 * t253, t232 * t238 - t236 * t250, t242 * t230 + t255 * t232 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end