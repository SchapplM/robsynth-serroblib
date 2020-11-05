% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRP5 (for one body)
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

function Tc_mdh = S6PRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t120 = cos(pkin(10));
	t119 = sin(pkin(10));
	t1 = [t120, -t119, 0, 0; t119, t120, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t121 = sin(pkin(10));
	t122 = sin(pkin(6));
	t130 = t121 * t122;
	t123 = cos(pkin(10));
	t129 = t123 * t122;
	t124 = cos(pkin(6));
	t125 = sin(qJ(2));
	t128 = t124 * t125;
	t126 = cos(qJ(2));
	t127 = t124 * t126;
	t1 = [-t121 * t128 + t123 * t126, -t121 * t127 - t123 * t125, t130, t123 * pkin(1) + pkin(7) * t130 + 0; t121 * t126 + t123 * t128, -t121 * t125 + t123 * t127, -t129, t121 * pkin(1) - pkin(7) * t129 + 0; t122 * t125, t122 * t126, t124, t124 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t134 = sin(pkin(6));
	t147 = t134 * pkin(7);
	t133 = sin(pkin(10));
	t136 = cos(pkin(6));
	t146 = t133 * t136;
	t137 = sin(qJ(3));
	t145 = t134 * t137;
	t139 = cos(qJ(3));
	t144 = t134 * t139;
	t135 = cos(pkin(10));
	t143 = t135 * t136;
	t138 = sin(qJ(2));
	t142 = t136 * t138;
	t140 = cos(qJ(2));
	t141 = t136 * t140;
	t132 = t133 * t140 + t135 * t142;
	t131 = t133 * t142 - t135 * t140;
	t1 = [-t131 * t139 + t133 * t145, t131 * t137 + t133 * t144, t133 * t141 + t135 * t138, (t135 * pkin(2) + pkin(8) * t146) * t140 + (-pkin(2) * t146 + t135 * pkin(8)) * t138 + t133 * t147 + t135 * pkin(1) + 0; t132 * t139 - t135 * t145, -t132 * t137 - t135 * t144, t133 * t138 - t135 * t141, (t133 * pkin(2) - pkin(8) * t143) * t140 + (pkin(2) * t143 + t133 * pkin(8)) * t138 - t135 * t147 + t133 * pkin(1) + 0; t136 * t137 + t138 * t144, t136 * t139 - t138 * t145, -t134 * t140, t136 * pkin(7) + qJ(1) + 0 + (pkin(2) * t138 - pkin(8) * t140) * t134; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->26), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->20)
	t151 = sin(pkin(6));
	t153 = cos(pkin(6));
	t155 = sin(qJ(2));
	t157 = cos(qJ(2));
	t154 = sin(qJ(3));
	t156 = cos(qJ(3));
	t161 = pkin(3) * t156 + qJ(4) * t154 + pkin(2);
	t159 = -pkin(8) * t157 + t161 * t155;
	t160 = t154 * pkin(3) - qJ(4) * t156 + pkin(7);
	t169 = t160 * t151 - t159 * t153;
	t165 = t153 * t155;
	t164 = t153 * t157;
	t163 = t154 * t151;
	t162 = t156 * t151;
	t158 = pkin(8) * t155 + t161 * t157 + pkin(1);
	t152 = cos(pkin(10));
	t150 = sin(pkin(10));
	t149 = t150 * t157 + t152 * t165;
	t148 = t150 * t165 - t152 * t157;
	t1 = [t150 * t164 + t152 * t155, t148 * t156 - t150 * t163, -t148 * t154 - t150 * t162, t169 * t150 + t158 * t152 + 0; t150 * t155 - t152 * t164, -t149 * t156 + t152 * t163, t149 * t154 + t152 * t162, t158 * t150 - t169 * t152 + 0; -t151 * t157, -t153 * t154 - t155 * t162, -t153 * t156 + t155 * t163, t159 * t151 + t160 * t153 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (70->34), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->31)
	t177 = sin(pkin(6));
	t179 = cos(pkin(6));
	t182 = sin(qJ(2));
	t185 = cos(qJ(2));
	t186 = pkin(4) + pkin(8);
	t181 = sin(qJ(3));
	t184 = cos(qJ(3));
	t187 = pkin(3) + pkin(9);
	t191 = qJ(4) * t181 + t187 * t184 + pkin(2);
	t189 = t191 * t182 - t186 * t185;
	t190 = qJ(4) * t184 - t187 * t181 - pkin(7);
	t202 = t190 * t177 + t189 * t179;
	t199 = t177 * t181;
	t198 = t177 * t185;
	t197 = t179 * t182;
	t196 = t179 * t185;
	t195 = t181 * t182;
	t194 = t181 * t185;
	t193 = t184 * t177;
	t188 = t186 * t182 + t191 * t185 + pkin(1);
	t183 = cos(qJ(5));
	t180 = sin(qJ(5));
	t178 = cos(pkin(10));
	t176 = sin(pkin(10));
	t175 = t177 * t195 - t179 * t184;
	t174 = t179 * t195 + t193;
	t173 = t176 * t196 + t178 * t182;
	t172 = t176 * t182 - t178 * t196;
	t171 = -t176 * t174 + t178 * t194;
	t170 = t178 * t174 + t176 * t194;
	t1 = [t171 * t180 + t173 * t183, t171 * t183 - t173 * t180, -(t176 * t197 - t178 * t185) * t184 + t176 * t199, -t202 * t176 + t188 * t178 + 0; t170 * t180 + t172 * t183, t170 * t183 - t172 * t180, (t176 * t185 + t178 * t197) * t184 - t178 * t199, t188 * t176 + t202 * t178 + 0; t175 * t180 - t183 * t198, t175 * t183 + t180 * t198, t179 * t181 + t182 * t193, t189 * t177 - t190 * t179 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:09:54
	% EndTime: 2020-11-04 21:09:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (96->38), mult. (174->60), div. (0->0), fcn. (208->10), ass. (0->32)
	t211 = sin(pkin(6));
	t213 = cos(pkin(6));
	t214 = sin(qJ(5));
	t217 = cos(qJ(5));
	t209 = t214 * pkin(5) - qJ(6) * t217 + qJ(4);
	t215 = sin(qJ(3));
	t218 = cos(qJ(3));
	t221 = pkin(3) + pkin(9);
	t224 = t209 * t218 - t221 * t215 - pkin(7);
	t216 = sin(qJ(2));
	t219 = cos(qJ(2));
	t225 = pkin(5) * t217 + qJ(6) * t214 + pkin(4) + pkin(8);
	t226 = t209 * t215 + t221 * t218 + pkin(2);
	t237 = t226 * t216 - t225 * t219;
	t238 = t224 * t211 + t237 * t213;
	t234 = t211 * t215;
	t233 = t211 * t219;
	t232 = t213 * t216;
	t231 = t213 * t219;
	t230 = t215 * t216;
	t229 = t215 * t219;
	t228 = t218 * t211;
	t222 = t225 * t216 + t226 * t219 + pkin(1);
	t212 = cos(pkin(10));
	t210 = sin(pkin(10));
	t208 = t211 * t230 - t213 * t218;
	t207 = t213 * t230 + t228;
	t206 = t210 * t231 + t212 * t216;
	t205 = t210 * t216 - t212 * t231;
	t204 = -t210 * t207 + t212 * t229;
	t203 = t212 * t207 + t210 * t229;
	t1 = [t204 * t214 + t206 * t217, -(t210 * t232 - t212 * t219) * t218 + t210 * t234, -t204 * t217 + t206 * t214, -t238 * t210 + t222 * t212 + 0; t203 * t214 + t205 * t217, (t210 * t219 + t212 * t232) * t218 - t212 * t234, -t203 * t217 + t205 * t214, t222 * t210 + t238 * t212 + 0; t208 * t214 - t217 * t233, t213 * t215 + t216 * t228, -t208 * t217 - t214 * t233, t237 * t211 - t224 * t213 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end