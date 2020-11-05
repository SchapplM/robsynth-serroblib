% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t103 = cos(pkin(9));
	t102 = sin(pkin(9));
	t1 = [t103, -t102, 0, 0; t102, t103, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t104 = sin(pkin(9));
	t105 = sin(pkin(5));
	t113 = t104 * t105;
	t106 = cos(pkin(9));
	t112 = t106 * t105;
	t107 = cos(pkin(5));
	t108 = sin(qJ(2));
	t111 = t107 * t108;
	t109 = cos(qJ(2));
	t110 = t107 * t109;
	t1 = [-t104 * t111 + t106 * t109, -t104 * t110 - t106 * t108, t113, t106 * pkin(1) + pkin(6) * t113 + 0; t104 * t109 + t106 * t111, -t104 * t108 + t106 * t110, -t112, t104 * pkin(1) - pkin(6) * t112 + 0; t105 * t108, t105 * t109, t107, t107 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t117 = sin(pkin(5));
	t130 = t117 * pkin(6);
	t116 = sin(pkin(9));
	t119 = cos(pkin(5));
	t129 = t116 * t119;
	t120 = sin(qJ(3));
	t128 = t117 * t120;
	t122 = cos(qJ(3));
	t127 = t117 * t122;
	t118 = cos(pkin(9));
	t126 = t118 * t119;
	t121 = sin(qJ(2));
	t125 = t119 * t121;
	t123 = cos(qJ(2));
	t124 = t119 * t123;
	t115 = -t116 * t125 + t118 * t123;
	t114 = t116 * t123 + t118 * t125;
	t1 = [t115 * t122 + t116 * t128, -t115 * t120 + t116 * t127, t116 * t124 + t118 * t121, (t118 * pkin(2) + pkin(7) * t129) * t123 + (-pkin(2) * t129 + t118 * pkin(7)) * t121 + t116 * t130 + t118 * pkin(1) + 0; t114 * t122 - t118 * t128, -t114 * t120 - t118 * t127, t116 * t121 - t118 * t124, (t116 * pkin(2) - pkin(7) * t126) * t123 + (pkin(2) * t126 + t116 * pkin(7)) * t121 - t118 * t130 + t116 * pkin(1) + 0; t119 * t120 + t121 * t127, t119 * t122 - t121 * t128, -t117 * t123, t119 * pkin(6) + qJ(1) + 0 + (pkin(2) * t121 - pkin(7) * t123) * t117; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t138 = sin(pkin(5));
	t140 = cos(pkin(5));
	t143 = sin(qJ(2));
	t146 = cos(qJ(2));
	t142 = sin(qJ(3));
	t145 = cos(qJ(3));
	t150 = t145 * pkin(3) + t142 * pkin(8) + pkin(2);
	t148 = -pkin(7) * t146 + t150 * t143;
	t149 = t142 * pkin(3) - t145 * pkin(8) + pkin(6);
	t161 = t149 * t138 - t148 * t140;
	t157 = t138 * t145;
	t156 = t138 * t146;
	t155 = t140 * t143;
	t154 = t140 * t146;
	t153 = t142 * t138;
	t152 = t143 * t145;
	t151 = t145 * t146;
	t147 = pkin(7) * t143 + t150 * t146 + pkin(1);
	t144 = cos(qJ(4));
	t141 = sin(qJ(4));
	t139 = cos(pkin(9));
	t137 = sin(pkin(9));
	t136 = t140 * t152 - t153;
	t135 = t138 * t152 + t140 * t142;
	t134 = t137 * t154 + t139 * t143;
	t133 = t137 * t143 - t139 * t154;
	t132 = -t137 * t136 + t139 * t151;
	t131 = t139 * t136 + t137 * t151;
	t1 = [t132 * t144 + t134 * t141, -t132 * t141 + t134 * t144, -(t137 * t155 - t139 * t146) * t142 - t137 * t157, t161 * t137 + t147 * t139 + 0; t131 * t144 + t133 * t141, -t131 * t141 + t133 * t144, (t137 * t146 + t139 * t155) * t142 + t139 * t157, t147 * t137 - t161 * t139 + 0; t135 * t144 - t141 * t156, -t135 * t141 - t144 * t156, -t140 * t145 + t143 * t153, t148 * t138 + t149 * t140 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:16
	% EndTime: 2020-11-04 20:07:16
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (83->36), mult. (180->60), div. (0->0), fcn. (214->10), ass. (0->31)
	t172 = sin(pkin(5));
	t174 = cos(pkin(5));
	t177 = sin(qJ(2));
	t180 = cos(qJ(2));
	t175 = sin(qJ(4));
	t178 = cos(qJ(4));
	t183 = t175 * pkin(4) - qJ(5) * t178 + pkin(7);
	t176 = sin(qJ(3));
	t179 = cos(qJ(3));
	t196 = pkin(4) * t178 + qJ(5) * t175 + pkin(3);
	t184 = t176 * pkin(8) + t179 * t196 + pkin(2);
	t182 = t184 * t177 - t183 * t180;
	t195 = -t179 * pkin(8) + t176 * t196 + pkin(6);
	t198 = t195 * t172 - t182 * t174;
	t193 = t172 * t179;
	t192 = t172 * t180;
	t191 = t174 * t177;
	t190 = t174 * t180;
	t189 = t176 * t172;
	t188 = t177 * t179;
	t187 = t179 * t180;
	t181 = t183 * t177 + t184 * t180 + pkin(1);
	t173 = cos(pkin(9));
	t171 = sin(pkin(9));
	t167 = t174 * t188 - t189;
	t166 = t172 * t188 + t174 * t176;
	t165 = t171 * t190 + t173 * t177;
	t164 = t171 * t177 - t173 * t190;
	t163 = -t171 * t167 + t173 * t187;
	t162 = t173 * t167 + t171 * t187;
	t1 = [t163 * t178 + t165 * t175, -(t171 * t191 - t173 * t180) * t176 - t171 * t193, t163 * t175 - t165 * t178, t198 * t171 + t181 * t173 + 0; t162 * t178 + t164 * t175, (t171 * t180 + t173 * t191) * t176 + t173 * t193, t162 * t175 - t164 * t178, t181 * t171 - t198 * t173 + 0; t166 * t178 - t175 * t192, -t174 * t179 + t177 * t189, t166 * t175 + t178 * t192, t182 * t172 + t195 * t174 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end