% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (116->8), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->16)
	t46 = pkin(9) + qJ(4);
	t44 = sin(t46);
	t45 = cos(t46);
	t47 = sin(pkin(7));
	t52 = cos(pkin(7)) * cos(pkin(8));
	t42 = t47 * t44 + t45 * t52;
	t39 = 0.1e1 / t42 ^ 2;
	t56 = qJD(4) * t39;
	t41 = t44 * t52 - t47 * t45;
	t38 = t41 ^ 2;
	t35 = t38 * t39 + 0.1e1;
	t53 = t42 * t56;
	t54 = t41 / t42 * t56;
	t55 = (t38 * t54 + t41 * t53) / t35 ^ 2;
	t33 = 0.1e1 / t35;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t55 + 0.2e1 * (t33 * t53 + (t33 * t54 - t39 * t55) * t41) * t41, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (1741->55), mult. (2271->132), div. (423->14), fcn. (2956->11), ass. (0->67)
	t142 = pkin(9) + qJ(4);
	t138 = sin(t142);
	t144 = sin(pkin(7));
	t145 = cos(pkin(8));
	t139 = cos(t142);
	t146 = cos(pkin(7));
	t166 = t139 * t146;
	t130 = t138 * t144 + t145 * t166;
	t147 = sin(qJ(5));
	t148 = cos(qJ(5));
	t143 = sin(pkin(8));
	t164 = t143 * t146;
	t155 = -t130 * t147 + t148 * t164;
	t174 = t155 * qJD(5);
	t162 = t146 * t138;
	t163 = t144 * t145;
	t128 = t139 * t163 - t162;
	t126 = t138 * t163 + t166;
	t165 = t143 * t138;
	t118 = atan2(-t126, t165);
	t114 = sin(t118);
	t115 = cos(t118);
	t101 = -t114 * t126 + t115 * t165;
	t98 = 0.1e1 / t101;
	t113 = t130 * t148 + t147 * t164;
	t109 = 0.1e1 / t113;
	t135 = 0.1e1 / t138;
	t99 = 0.1e1 / t101 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t136 = 0.1e1 / t138 ^ 2;
	t108 = t155 ^ 2;
	t105 = t108 * t110 + 0.1e1;
	t129 = -t144 * t139 + t145 * t162;
	t122 = t129 * qJD(4);
	t106 = t113 * qJD(5) - t122 * t147;
	t169 = t110 * t155;
	t107 = -t122 * t148 + t174;
	t170 = t107 * t109 * t110;
	t173 = 0.1e1 / t105 ^ 2 * (-t106 * t169 - t108 * t170);
	t167 = t136 * t139;
	t159 = t126 * t167;
	t156 = -t128 * t135 + t159;
	t124 = t126 ^ 2;
	t141 = 0.1e1 / t143 ^ 2;
	t119 = t124 * t136 * t141 + 0.1e1;
	t116 = 0.1e1 / t119;
	t140 = 0.1e1 / t143;
	t168 = t116 * t140;
	t94 = t156 * t168;
	t172 = t126 * t94;
	t171 = t129 * t99;
	t161 = qJD(4) * t139;
	t160 = -0.2e1 * t173;
	t157 = -t109 * t147 - t148 * t169;
	t137 = t135 * t136;
	t125 = t129 ^ 2;
	t123 = t130 * qJD(4);
	t121 = t128 * qJD(4);
	t120 = t126 * qJD(4);
	t103 = 0.1e1 / t105;
	t100 = t98 * t99;
	t97 = t125 * t99 + 0.1e1;
	t93 = (qJD(4) * t159 - t121 * t135) * t168;
	t91 = (t139 * t143 - t172) * t115 + (-t94 * t165 - t128) * t114;
	t90 = (-t126 * t93 + t143 * t161) * t115 + (-t93 * t165 - t121) * t114;
	t89 = (-0.2e1 * t156 / t119 ^ 2 * (t121 * t126 * t136 - t124 * t137 * t161) * t141 + (t121 * t167 + t120 * t135 + (t128 * t167 + (-0.2e1 * t137 * t139 ^ 2 - t135) * t126) * qJD(4)) * t116) * t140;
	t1 = [0, 0, 0, t89, 0; 0, 0, 0, 0.2e1 * (-t130 * t98 + t91 * t171) / t97 ^ 2 * (-t100 * t125 * t90 + t123 * t171) + (-t91 * t123 * t99 - t122 * t98 + (0.2e1 * t100 * t129 * t91 - t130 * t99) * t90 + (-(-t121 * t94 - t126 * t89 - t128 * t93 + (-t93 * t94 - qJD(4)) * t165) * t115 - (t93 * t172 + t120 + (-t138 * t89 + (-qJD(4) * t94 - t93) * t139) * t143) * t114) * t171) / t97, 0; 0, 0, 0, t157 * t129 * t160 + (t157 * t123 + ((-qJD(5) * t109 + 0.2e1 * t155 * t170) * t148 + (t106 * t148 + (t107 + t174) * t147) * t110) * t129) * t103, t160 - 0.2e1 * (t103 * t106 * t110 - (-t103 * t170 - t110 * t173) * t155) * t155;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end