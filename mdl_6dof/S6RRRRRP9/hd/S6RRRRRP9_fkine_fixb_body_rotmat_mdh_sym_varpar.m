% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t103 = cos(qJ(1));
	t102 = sin(qJ(1));
	t1 = [t103, -t102, 0, 0; t102, t103, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t104 = sin(pkin(6));
	t107 = sin(qJ(1));
	t115 = t107 * t104;
	t106 = sin(qJ(2));
	t114 = t107 * t106;
	t108 = cos(qJ(2));
	t113 = t107 * t108;
	t109 = cos(qJ(1));
	t112 = t109 * t104;
	t111 = t109 * t106;
	t110 = t109 * t108;
	t105 = cos(pkin(6));
	t1 = [-t105 * t114 + t110, -t105 * t113 - t111, t115, t109 * pkin(1) + pkin(8) * t115 + 0; t105 * t111 + t113, t105 * t110 - t114, -t112, t107 * pkin(1) - pkin(8) * t112 + 0; t104 * t106, t104 * t108, t105, t105 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (29->23), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->20)
	t119 = sin(pkin(6));
	t124 = cos(qJ(3));
	t134 = t119 * t124;
	t121 = sin(qJ(3));
	t133 = t121 * t119;
	t122 = sin(qJ(2));
	t132 = t122 * t124;
	t123 = sin(qJ(1));
	t131 = t123 * t122;
	t125 = cos(qJ(2));
	t130 = t123 * t125;
	t126 = cos(qJ(1));
	t129 = t126 * t122;
	t128 = t126 * t125;
	t127 = pkin(2) * t122 - pkin(9) * t125;
	t120 = cos(pkin(6));
	t118 = t125 * pkin(2) + t122 * pkin(9) + pkin(1);
	t117 = -t120 * t129 - t130;
	t116 = t119 * pkin(8) - t127 * t120;
	t1 = [(-t120 * t132 + t133) * t123 + t124 * t128, (t120 * t131 - t128) * t121 + t123 * t134, t120 * t130 + t129, t116 * t123 + t118 * t126 + 0; -t117 * t124 - t126 * t133, t117 * t121 - t126 * t134, -t120 * t128 + t131, -t116 * t126 + t118 * t123 + 0; t119 * t132 + t120 * t121, t120 * t124 - t122 * t133, -t119 * t125, t120 * pkin(8) + t127 * t119 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->34), mult. (112->60), div. (0->0), fcn. (146->10), ass. (0->29)
	t142 = sin(pkin(6));
	t149 = cos(qJ(3));
	t163 = t142 * t149;
	t144 = sin(qJ(4));
	t150 = cos(qJ(2));
	t162 = t144 * t150;
	t145 = sin(qJ(3));
	t161 = t145 * t142;
	t146 = sin(qJ(2));
	t160 = t146 * t149;
	t147 = sin(qJ(1));
	t159 = t147 * t150;
	t148 = cos(qJ(4));
	t158 = t148 * t150;
	t157 = t149 * t150;
	t151 = cos(qJ(1));
	t156 = t151 * t146;
	t155 = t151 * t150;
	t140 = t149 * pkin(3) + t145 * pkin(10) + pkin(2);
	t154 = t150 * pkin(9) - t140 * t146;
	t153 = t145 * pkin(3) - t149 * pkin(10) + pkin(8);
	t143 = cos(pkin(6));
	t138 = t143 * t160 - t161;
	t152 = t138 * t144 + t143 * t158;
	t139 = t144 * t157 - t148 * t146;
	t137 = t142 * t160 + t143 * t145;
	t136 = t146 * pkin(9) + t140 * t150 + pkin(1);
	t135 = t142 * t153 + t154 * t143;
	t1 = [(-t138 * t147 + t149 * t155) * t148 + (t143 * t159 + t156) * t144, -t151 * t139 + t152 * t147, (-t147 * t143 * t146 + t155) * t145 - t147 * t163, t135 * t147 + t136 * t151 + 0; (t138 * t148 - t143 * t162) * t151 + t147 * (t146 * t144 + t148 * t157), -t147 * t139 - t152 * t151, (t143 * t156 + t159) * t145 + t151 * t163, -t135 * t151 + t136 * t147 + 0; t137 * t148 - t142 * t162, -t137 * t144 - t142 * t158, -t143 * t149 + t146 * t161, -t154 * t142 + t153 * t143 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (90->36), mult. (125->55), div. (0->0), fcn. (159->12), ass. (0->34)
	t177 = sin(pkin(6));
	t182 = cos(qJ(3));
	t197 = t177 * t182;
	t183 = cos(qJ(2));
	t196 = t177 * t183;
	t179 = sin(qJ(3));
	t195 = t179 * t177;
	t180 = sin(qJ(2));
	t194 = t180 * t182;
	t181 = sin(qJ(1));
	t193 = t181 * t180;
	t192 = t181 * t183;
	t184 = cos(qJ(1));
	t191 = t184 * t180;
	t190 = t184 * t183;
	t173 = cos(qJ(4)) * pkin(4) + pkin(3);
	t185 = pkin(11) + pkin(10);
	t166 = t173 * t182 + t185 * t179 + pkin(2);
	t172 = sin(qJ(4)) * pkin(4) + pkin(9);
	t189 = t166 * t180 - t172 * t183;
	t188 = t173 * t179 - t185 * t182 + pkin(8);
	t178 = cos(pkin(6));
	t168 = t178 * t194 - t195;
	t187 = -t168 * t181 + t182 * t190;
	t186 = t168 * t184 + t182 * t192;
	t176 = qJ(4) + qJ(5);
	t175 = cos(t176);
	t174 = sin(t176);
	t170 = t178 * t192 + t191;
	t169 = t178 * t190 - t193;
	t167 = t177 * t194 + t178 * t179;
	t165 = t166 * t183 + t172 * t180 + pkin(1);
	t164 = t188 * t177 - t189 * t178;
	t1 = [t170 * t174 + t187 * t175, t170 * t175 - t187 * t174, (-t178 * t193 + t190) * t179 - t181 * t197, t164 * t181 + t165 * t184 + 0; -t174 * t169 + t186 * t175, -t175 * t169 - t186 * t174, (t178 * t191 + t192) * t179 + t184 * t197, -t164 * t184 + t165 * t181 + 0; t167 * t175 - t174 * t196, -t167 * t174 - t175 * t196, -t178 * t182 + t180 * t195, t189 * t177 + t188 * t178 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:44:56
	% EndTime: 2020-11-04 22:44:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (95->45), mult. (152->62), div. (0->0), fcn. (198->12), ass. (0->35)
	t218 = sin(qJ(2));
	t231 = pkin(2) * t218;
	t215 = sin(pkin(6));
	t220 = cos(qJ(3));
	t230 = t215 * t220;
	t221 = cos(qJ(2));
	t229 = t215 * t221;
	t217 = sin(qJ(3));
	t228 = t217 * t215;
	t227 = t218 * t220;
	t219 = sin(qJ(1));
	t226 = t219 * t218;
	t225 = t219 * t221;
	t222 = cos(qJ(1));
	t224 = t222 * t218;
	t223 = t222 * t221;
	t216 = cos(pkin(6));
	t214 = qJ(4) + qJ(5);
	t213 = -qJ(6) - pkin(11) - pkin(10);
	t212 = cos(t214);
	t211 = sin(t214);
	t210 = t221 * pkin(2) + t218 * pkin(9) + pkin(1);
	t209 = pkin(5) * t211 + sin(qJ(4)) * pkin(4);
	t208 = pkin(5) * t212 + cos(qJ(4)) * pkin(4) + pkin(3);
	t207 = t216 * t225 + t224;
	t206 = -t216 * t223 + t226;
	t205 = -t216 * t224 - t225;
	t204 = t215 * t227 + t216 * t217;
	t203 = -t216 * t220 + t218 * t228;
	t202 = t215 * pkin(8) + (pkin(9) * t221 - t231) * t216;
	t201 = -t205 * t220 - t222 * t228;
	t200 = t205 * t217 - t222 * t230;
	t199 = (-t216 * t227 + t228) * t219 + t220 * t223;
	t198 = (t216 * t226 - t223) * t217 + t219 * t230;
	t1 = [t199 * t212 + t207 * t211, -t199 * t211 + t207 * t212, -t198, t198 * t213 + t199 * t208 + t202 * t219 + t207 * t209 + t210 * t222 + 0; t201 * t212 + t206 * t211, -t201 * t211 + t206 * t212, -t200, t200 * t213 + t201 * t208 - t202 * t222 + t206 * t209 + t210 * t219 + 0; t204 * t212 - t211 * t229, -t204 * t211 - t212 * t229, t203, t216 * pkin(8) - t203 * t213 + t204 * t208 + pkin(7) + 0 + (t231 + (-pkin(9) - t209) * t221) * t215; 0, 0, 0, 1;];
	Tc_mdh = t1;
end