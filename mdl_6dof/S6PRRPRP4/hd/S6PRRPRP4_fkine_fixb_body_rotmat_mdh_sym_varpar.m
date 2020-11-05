% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t113 = cos(pkin(10));
	t112 = sin(pkin(10));
	t1 = [t113, -t112, 0, 0; t112, t113, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t114 = sin(pkin(10));
	t115 = sin(pkin(6));
	t123 = t114 * t115;
	t116 = cos(pkin(10));
	t122 = t116 * t115;
	t117 = cos(pkin(6));
	t118 = sin(qJ(2));
	t121 = t117 * t118;
	t119 = cos(qJ(2));
	t120 = t117 * t119;
	t1 = [-t114 * t121 + t116 * t119, -t114 * t120 - t116 * t118, t123, t116 * pkin(1) + pkin(7) * t123 + 0; t114 * t119 + t116 * t121, -t114 * t118 + t116 * t120, -t122, t114 * pkin(1) - pkin(7) * t122 + 0; t115 * t118, t115 * t119, t117, t117 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t127 = sin(pkin(6));
	t140 = t127 * pkin(7);
	t126 = sin(pkin(10));
	t129 = cos(pkin(6));
	t139 = t126 * t129;
	t130 = sin(qJ(3));
	t138 = t127 * t130;
	t132 = cos(qJ(3));
	t137 = t127 * t132;
	t128 = cos(pkin(10));
	t136 = t128 * t129;
	t131 = sin(qJ(2));
	t135 = t129 * t131;
	t133 = cos(qJ(2));
	t134 = t129 * t133;
	t125 = t126 * t133 + t128 * t135;
	t124 = t126 * t135 - t128 * t133;
	t1 = [-t124 * t132 + t126 * t138, t124 * t130 + t126 * t137, t126 * t134 + t128 * t131, (t128 * pkin(2) + pkin(8) * t139) * t133 + (-pkin(2) * t139 + t128 * pkin(8)) * t131 + t126 * t140 + t128 * pkin(1) + 0; t125 * t132 - t128 * t138, -t125 * t130 - t128 * t137, t126 * t131 - t128 * t134, (t126 * pkin(2) - pkin(8) * t136) * t133 + (pkin(2) * t136 + t126 * pkin(8)) * t131 - t128 * t140 + t126 * pkin(1) + 0; t129 * t130 + t131 * t137, t129 * t132 - t131 * t138, -t127 * t133, t129 * pkin(7) + qJ(1) + 0 + (pkin(2) * t131 - pkin(8) * t133) * t127; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t144 = sin(pkin(6));
	t146 = cos(pkin(6));
	t148 = sin(qJ(2));
	t150 = cos(qJ(2));
	t147 = sin(qJ(3));
	t149 = cos(qJ(3));
	t154 = pkin(3) * t149 + qJ(4) * t147 + pkin(2);
	t152 = -pkin(8) * t150 + t154 * t148;
	t153 = pkin(3) * t147 - qJ(4) * t149 + pkin(7);
	t162 = t153 * t144 - t152 * t146;
	t158 = t144 * t147;
	t157 = t144 * t149;
	t156 = t146 * t148;
	t155 = t146 * t150;
	t151 = pkin(8) * t148 + t154 * t150 + pkin(1);
	t145 = cos(pkin(10));
	t143 = sin(pkin(10));
	t142 = t143 * t150 + t145 * t156;
	t141 = t143 * t156 - t145 * t150;
	t1 = [t143 * t155 + t145 * t148, t141 * t149 - t143 * t158, -t141 * t147 - t143 * t157, t162 * t143 + t151 * t145 + 0; t143 * t148 - t145 * t155, -t142 * t149 + t145 * t158, t142 * t147 + t145 * t157, t151 * t143 - t162 * t145 + 0; -t144 * t150, -t146 * t147 - t148 * t157, -t146 * t149 + t148 * t158, t152 * t144 + t153 * t146 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (70->34), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->31)
	t170 = sin(pkin(6));
	t172 = cos(pkin(6));
	t175 = sin(qJ(2));
	t178 = cos(qJ(2));
	t179 = pkin(4) + pkin(8);
	t174 = sin(qJ(3));
	t177 = cos(qJ(3));
	t180 = pkin(3) + pkin(9);
	t184 = qJ(4) * t174 + t180 * t177 + pkin(2);
	t182 = t184 * t175 - t179 * t178;
	t183 = qJ(4) * t177 - t180 * t174 - pkin(7);
	t195 = t183 * t170 + t182 * t172;
	t192 = t170 * t174;
	t191 = t170 * t178;
	t190 = t172 * t175;
	t189 = t172 * t178;
	t188 = t174 * t175;
	t187 = t174 * t178;
	t186 = t177 * t170;
	t181 = t179 * t175 + t184 * t178 + pkin(1);
	t176 = cos(qJ(5));
	t173 = sin(qJ(5));
	t171 = cos(pkin(10));
	t169 = sin(pkin(10));
	t168 = t170 * t188 - t172 * t177;
	t167 = t172 * t188 + t186;
	t166 = t169 * t189 + t171 * t175;
	t165 = t169 * t175 - t171 * t189;
	t164 = -t169 * t167 + t171 * t187;
	t163 = t171 * t167 + t169 * t187;
	t1 = [t164 * t173 + t166 * t176, t164 * t176 - t166 * t173, -(t169 * t190 - t171 * t178) * t177 + t169 * t192, -t195 * t169 + t181 * t171 + 0; t163 * t173 + t165 * t176, t163 * t176 - t165 * t173, (t169 * t178 + t171 * t190) * t177 - t171 * t192, t181 * t169 + t195 * t171 + 0; t168 * t173 - t176 * t191, t168 * t176 + t173 * t191, t172 * t174 + t175 * t186, t182 * t170 - t183 * t172 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:31
	% EndTime: 2020-11-04 21:09:31
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (91->36), mult. (155->58), div. (0->0), fcn. (189->10), ass. (0->32)
	t205 = sin(pkin(6));
	t207 = cos(pkin(6));
	t208 = sin(qJ(5));
	t202 = t208 * pkin(5) + qJ(4);
	t203 = qJ(6) + pkin(3) + pkin(9);
	t209 = sin(qJ(3));
	t212 = cos(qJ(3));
	t217 = t202 * t212 - t203 * t209 - pkin(7);
	t210 = sin(qJ(2));
	t213 = cos(qJ(2));
	t218 = t202 * t209 + t203 * t212 + pkin(2);
	t211 = cos(qJ(5));
	t219 = t211 * pkin(5) + pkin(4) + pkin(8);
	t230 = -t218 * t210 + t219 * t213;
	t231 = -t217 * t205 + t230 * t207;
	t226 = t205 * t209;
	t225 = t205 * t213;
	t224 = t207 * t210;
	t223 = t207 * t213;
	t222 = t209 * t210;
	t221 = t209 * t213;
	t220 = t212 * t205;
	t215 = t219 * t210 + t218 * t213 + pkin(1);
	t206 = cos(pkin(10));
	t204 = sin(pkin(10));
	t201 = t205 * t222 - t207 * t212;
	t200 = t207 * t222 + t220;
	t199 = t204 * t223 + t206 * t210;
	t198 = t204 * t210 - t206 * t223;
	t197 = -t204 * t200 + t206 * t221;
	t196 = t206 * t200 + t204 * t221;
	t1 = [t197 * t208 + t199 * t211, t197 * t211 - t199 * t208, -(t204 * t224 - t206 * t213) * t212 + t204 * t226, t231 * t204 + t215 * t206 + 0; t196 * t208 + t198 * t211, t196 * t211 - t198 * t208, (t204 * t213 + t206 * t224) * t212 - t206 * t226, t215 * t204 - t231 * t206 + 0; t201 * t208 - t211 * t225, t201 * t211 + t208 * t225, t207 * t209 + t210 * t220, -t230 * t205 - t217 * t207 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end