% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR5
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
%   Wie in S5PRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:30
	% EndTime: 2019-12-29 15:31:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:20
	% EndTime: 2019-12-29 15:31:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:24
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (309->29), mult. (767->69), div. (195->13), fcn. (861->7), ass. (0->42)
	t63 = sin(qJ(2));
	t57 = t63 ^ 2;
	t64 = cos(qJ(2));
	t59 = 0.1e1 / t64 ^ 2;
	t82 = t57 * t59;
	t58 = 0.1e1 / t64;
	t61 = sin(pkin(7));
	t53 = t61 ^ 2;
	t90 = t53 * t59;
	t89 = t58 * t82;
	t88 = qJD(2) * t61;
	t80 = t61 * t63;
	t47 = atan2(-t80, -t64);
	t45 = sin(t47);
	t46 = cos(t47);
	t41 = -t45 * t80 - t46 * t64;
	t38 = 0.1e1 / t41;
	t39 = 0.1e1 / t41 ^ 2;
	t52 = t53 * t82 + 0.1e1;
	t49 = 0.1e1 / t52;
	t76 = 0.1e1 + t82;
	t43 = t76 * t61 * t49;
	t83 = t45 * t64;
	t74 = t46 * t63 - t61 * t83;
	t78 = t46 * t80;
	t75 = -t78 + t83;
	t33 = t75 * t43 + t74;
	t87 = 0.2e1 * t33;
	t62 = cos(pkin(7));
	t54 = t62 ^ 2;
	t37 = t39 * t54 * t57 + 0.1e1;
	t84 = t39 * t64;
	t42 = qJD(2) * t43;
	t32 = t74 * qJD(2) + t75 * t42;
	t85 = t32 * t38 * t39;
	t86 = (qJD(2) * t63 * t84 - t57 * t85) * t54 / t37 ^ 2;
	t79 = t43 - t61;
	t77 = t43 * t61 - 0.1e1;
	t51 = 0.1e1 + 0.1e1 / t62 ^ 2 * t90;
	t35 = 0.1e1 / t37;
	t34 = 0.2e1 * (t49 - t76 / t52 ^ 2 * t53) * (t58 + t89) * t63 * t88;
	t1 = [0, t34, 0, 0, 0; 0, ((-0.2e1 * t38 * t86 + (-qJD(2) * t33 - t32) * t39 * t35) * t64 + (t39 * t86 * t87 + (t85 * t87 - qJD(2) * t38 - (t34 * t45 + (-t77 * qJD(2) + t79 * t42) * t46) * t84 + (t34 * t78 - (-t79 * qJD(2) + t77 * t42) * t63 * t45) * t39) * t35) * t63) * t62, 0, 0, 0; 0, (0.2e1 / t54 / t51 ^ 2 * t89 * t90 + (-t58 - 0.2e1 * t89) / t51) / t62 * t88, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:26
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (981->53), mult. (3022->121), div. (281->12), fcn. (3908->11), ass. (0->63)
	t155 = sin(pkin(8));
	t159 = cos(qJ(2));
	t184 = cos(pkin(8));
	t185 = sin(qJ(2));
	t148 = t185 * t155 + t159 * t184;
	t156 = cos(pkin(7));
	t141 = t148 * t156;
	t157 = sin(qJ(5));
	t158 = cos(qJ(5));
	t183 = sin(pkin(7));
	t167 = -t141 * t157 - t183 * t158;
	t187 = t167 * qJD(5);
	t171 = t184 * t183;
	t173 = t155 * t183;
	t139 = t159 * t173 - t185 * t171;
	t125 = atan2(-t139, t148);
	t120 = sin(t125);
	t121 = cos(t125);
	t113 = -t120 * t139 + t121 * t148;
	t110 = 0.1e1 / t113;
	t131 = t141 * t158 - t183 * t157;
	t127 = 0.1e1 / t131;
	t145 = 0.1e1 / t148;
	t111 = 0.1e1 / t113 ^ 2;
	t128 = 0.1e1 / t131 ^ 2;
	t146 = 0.1e1 / t148 ^ 2;
	t149 = t159 * t155 - t185 * t184;
	t142 = t149 * t156;
	t186 = 0.2e1 * t142;
	t136 = t139 ^ 2;
	t124 = t136 * t146 + 0.1e1;
	t122 = 0.1e1 / t124;
	t166 = t159 * t171 + t185 * t173;
	t132 = t166 * qJD(2);
	t144 = t149 * qJD(2);
	t177 = t139 * t146;
	t105 = (t132 * t145 + t144 * t177) * t122;
	t170 = -t120 * t148 - t121 * t139;
	t102 = t170 * t105 + t120 * t132 + t121 * t144;
	t182 = t102 * t110 * t111;
	t126 = t167 ^ 2;
	t117 = t126 * t128 + 0.1e1;
	t135 = t156 * t144;
	t118 = t131 * qJD(5) + t135 * t157;
	t178 = t128 * t167;
	t119 = t135 * t158 + t187;
	t179 = t119 * t127 * t128;
	t181 = (-t118 * t178 - t126 * t179) / t117 ^ 2;
	t180 = t111 * t142;
	t176 = t139 * t149;
	t175 = t144 * t145 * t146;
	t169 = -t127 * t157 - t158 * t178;
	t168 = t145 * t166 + t146 * t176;
	t143 = t148 * qJD(2);
	t137 = t142 ^ 2;
	t134 = t156 * t143;
	t133 = t139 * qJD(2);
	t115 = 0.1e1 / t117;
	t109 = t137 * t111 + 0.1e1;
	t106 = t168 * t122;
	t103 = t170 * t106 + t120 * t166 + t121 * t149;
	t101 = -0.2e1 * t168 / t124 ^ 2 * (-t132 * t177 - t136 * t175) + (-0.2e1 * t175 * t176 + t133 * t145 + (-t132 * t149 - t139 * t143 - t144 * t166) * t146) * t122;
	t1 = [0, t101, 0, 0, 0; 0, 0.2e1 * (t103 * t180 + t110 * t141) / t109 ^ 2 * (-t134 * t180 - t137 * t182) + (t103 * t182 * t186 - t135 * t110 + (t141 * t102 + t103 * t134 + (-(-t101 * t139 + t106 * t132 - t143 + (-t106 * t148 + t166) * t105) * t121 - (-t101 * t148 - t106 * t144 + t133 + (t106 * t139 - t149) * t105) * t120) * t142) * t111) / t109, 0, 0, 0; 0, t169 * t181 * t186 + (t169 * t134 + ((qJD(5) * t127 - 0.2e1 * t167 * t179) * t158 + (-t118 * t158 + (-t119 - t187) * t157) * t128) * t142) * t115, 0, 0, -0.2e1 * t181 - 0.2e1 * (t115 * t118 * t128 - (-t115 * t179 - t128 * t181) * t167) * t167;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end