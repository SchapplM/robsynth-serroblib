% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP5
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
%   Wie in S5RPPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:07
	% EndTime: 2019-12-29 16:05:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:13
	% EndTime: 2019-12-29 16:05:13
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (213->24), mult. (456->71), div. (110->13), fcn. (596->7), ass. (0->33)
	t68 = sin(qJ(1));
	t80 = qJD(1) * t68;
	t67 = cos(pkin(7));
	t58 = 0.1e1 / t67 ^ 2;
	t66 = sin(pkin(7));
	t56 = t66 ^ 2;
	t61 = t68 ^ 2;
	t53 = t56 * t58 * t61 + 0.1e1;
	t69 = cos(qJ(1));
	t62 = t69 ^ 2;
	t84 = 0.1e1 / t53 ^ 2 * t62;
	t89 = t58 * t84;
	t49 = 0.1e1 / t53;
	t87 = (t49 - 0.1e1) * t66;
	t81 = t68 * t66;
	t48 = atan2(-t81, -t67);
	t46 = sin(t48);
	t47 = cos(t48);
	t44 = -t46 * t81 - t47 * t67;
	t41 = 0.1e1 / t44;
	t57 = 0.1e1 / t67;
	t42 = 0.1e1 / t44 ^ 2;
	t85 = t42 * t69;
	t83 = t56 * t57;
	t82 = t61 / t69 ^ 2;
	t79 = t57 * t89;
	t37 = (-t47 * t49 * t68 * t83 + t46 * t87) * t69;
	t55 = t66 * t56;
	t54 = t58 * t82 + 0.1e1;
	t43 = t41 * t42;
	t40 = t42 * t56 * t62 + 0.1e1;
	t36 = qJD(1) * t37;
	t1 = [(-t49 * t57 * t66 - 0.2e1 * t55 * t79) * t80, 0, 0, 0, 0; (0.2e1 * (t37 * t85 + t41 * t68) / t40 ^ 2 * (-t36 * t43 * t62 - t80 * t85) * t56 + ((0.2e1 * t37 * t43 * t69 + t42 * t68) * t36 + (-t69 * t41 + ((t37 + (t55 * t89 + t87) * t69 * t46) * t68 - (0.2e1 * t61 * t56 ^ 2 * t79 + (t84 + (t61 - 0.2e1 * t62) * t49) * t83) * t69 * t47) * t42) * qJD(1)) / t40) * t66, 0, 0, 0, 0; 0.2e1 * (-0.1e1 / t54 + (0.1e1 + t82) / t54 ^ 2 * t58) * (t61 / t62 + 0.1e1) / t69 * t57 * t80, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:13
	% EndTime: 2019-12-29 16:05:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:14
	% DurationCPUTime: 1.34s
	% Computational Cost: add. (1478->65), mult. (4561->146), div. (413->12), fcn. (5923->9), ass. (0->67)
	t124 = cos(pkin(7));
	t125 = sin(qJ(4));
	t159 = sin(pkin(7));
	t168 = cos(qJ(4));
	t115 = -t124 * t125 + t159 * t168;
	t114 = t124 * t168 + t159 * t125;
	t111 = 0.1e1 / t114;
	t112 = 0.1e1 / t114 ^ 2;
	t110 = t115 * qJD(4);
	t155 = t110 * t112;
	t154 = t111 * t155;
	t126 = sin(qJ(1));
	t171 = t115 * t126;
	t106 = t114 * t126;
	t169 = cos(qJ(1));
	t133 = t115 * t169;
	t88 = -qJD(1) * t133 + qJD(4) * t106;
	t99 = t171 ^ 2;
	t94 = t112 * t99 + 0.1e1;
	t166 = (-t112 * t171 * t88 - t99 * t154) / t94 ^ 2;
	t172 = 0.2e1 * t111 * t166;
	t108 = t114 * t169;
	t157 = t171 * t111;
	t95 = atan2(t171, t114);
	t90 = sin(t95);
	t91 = cos(t95);
	t92 = 0.1e1 / t94;
	t170 = (-t91 * t157 + t90) * t92 - t90;
	t84 = t114 * t91 + t171 * t90;
	t81 = 0.1e1 / t84;
	t101 = 0.1e1 / t108;
	t102 = 0.1e1 / t108 ^ 2;
	t82 = 0.1e1 / t84 ^ 2;
	t140 = -t114 * t90 + t171 * t91;
	t136 = -t111 * t88 - t155 * t171;
	t75 = t136 * t92;
	t72 = t110 * t91 + t140 * t75 - t88 * t90;
	t167 = t72 * t81 * t82;
	t123 = t126 ^ 2;
	t145 = qJD(1) * t169;
	t158 = t102 * t126;
	t87 = -qJD(1) * t106 + t133 * qJD(4);
	t164 = t101 * t102 * t87;
	t98 = t102 * t123 + 0.1e1;
	t165 = (-t123 * t164 + t145 * t158) / t98 ^ 2;
	t163 = t133 * t82;
	t162 = t133 * t90;
	t161 = t133 * t91;
	t96 = 0.1e1 / t98;
	t160 = t126 * t96;
	t156 = t171 * t115;
	t100 = t133 ^ 2;
	t79 = t100 * t82 + 0.1e1;
	t86 = -qJD(1) * t171 - qJD(4) * t108;
	t152 = 0.2e1 * (-t100 * t167 + t86 * t163) / t79 ^ 2;
	t151 = 0.2e1 * t167;
	t150 = -0.2e1 * t166;
	t149 = -0.2e1 * t164;
	t135 = -t106 * t111 - t112 * t156;
	t109 = t114 * qJD(4);
	t89 = t108 * qJD(1) + qJD(4) * t171;
	t77 = 0.1e1 / t79;
	t76 = t135 * t92;
	t74 = t170 * t133;
	t73 = -t106 * t90 + t115 * t91 + t140 * t76;
	t71 = t135 * t150 + (0.2e1 * t154 * t156 - t111 * t89 + (t106 * t110 + t109 * t171 + t115 * t88) * t112) * t92;
	t1 = [t111 * t86 * t92 - (t92 * t155 + t172) * t133, 0, 0, t71, 0; -t171 * t81 * t152 + (-t88 * t81 + (-t171 * t72 - t74 * t86) * t82) * t77 - (-t74 * t151 * t77 + (-t74 * t152 + ((t75 * t92 * t157 + t150) * t162 + (t171 * t172 - t75 + (-t136 + t75) * t92) * t161 + t170 * t86) * t77) * t82) * t133, 0, 0, (-t108 * t81 - t73 * t163) * t152 + (-t73 * t133 * t151 + t87 * t81 + (-t108 * t72 + t73 * t86 + (t171 * t71 - t76 * t88 - t109 + (-t114 * t76 - t106) * t75) * t161 + (-t110 * t76 - t114 * t71 - t89 + (-t171 * t76 - t115) * t75) * t162) * t82) * t77, 0; 0.2e1 * (-t169 * t101 - t106 * t158) * t165 + ((-qJD(1) * t101 + t106 * t149) * t126 + (t106 * t145 + t126 * t89 - t169 * t87) * t102) * t96, 0, 0, -t133 * t149 * t160 + (-t86 * t160 - (-0.2e1 * t126 * t165 + t96 * t145) * t133) * t102, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end