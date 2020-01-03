% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP5
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
%   Wie in S5RPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (118->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:21
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (3227->71), mult. (2645->160), div. (674->13), fcn. (3179->7), ass. (0->78)
	t104 = qJ(1) + pkin(8) + qJ(3);
	t102 = sin(t104);
	t103 = cos(t104);
	t99 = 0.1e1 / t103;
	t151 = t102 * t99;
	t139 = qJD(4) * t102;
	t112 = cos(qJ(4));
	t108 = 0.1e1 / t112;
	t111 = sin(qJ(4));
	t107 = t111 ^ 2;
	t109 = 0.1e1 / t112 ^ 2;
	t142 = t107 * t109;
	t129 = t102 * t142;
	t138 = qJD(4) * t112;
	t105 = qJD(1) + qJD(3);
	t143 = t105 * t111;
	t97 = t102 ^ 2;
	t95 = t97 * t142 + 0.1e1;
	t93 = 0.1e1 / t95;
	t76 = (-(-t102 * t138 - t103 * t143) * t108 + qJD(4) * t129) * t93;
	t163 = t76 - t139;
	t98 = t103 ^ 2;
	t162 = t105 * (0.1e1 / t98 * t97 + 0.1e1) * t151;
	t145 = t102 * t111;
	t92 = atan2(-t145, -t112);
	t90 = sin(t92);
	t132 = t90 * t145;
	t91 = cos(t92);
	t85 = -t112 * t91 - t132;
	t82 = 0.1e1 / t85;
	t83 = 0.1e1 / t85 ^ 2;
	t160 = 0.2e1 * t111;
	t159 = t93 - 0.1e1;
	t130 = t83 * t138;
	t146 = t102 * t105;
	t134 = t83 * t146;
	t144 = t103 * t105;
	t149 = t112 * t90;
	t152 = t102 * t76;
	t71 = t163 * t149 + (-t90 * t144 + (qJD(4) - t152) * t91) * t111;
	t157 = t71 * t82 * t83;
	t80 = t107 * t83 * t98 + 0.1e1;
	t158 = (t98 * t111 * t130 + (-t103 * t134 - t98 * t157) * t107) / t80 ^ 2;
	t78 = 0.1e1 / t80;
	t156 = t78 * t83;
	t106 = t111 * t107;
	t110 = t108 * t109;
	t141 = t108 * t111;
	t120 = qJD(4) * (t106 * t110 + t141);
	t155 = (t97 * t120 + t129 * t144) / t95 ^ 2;
	t154 = t82 * t78;
	t153 = 0.1e1 / t103 ^ 2 * t97;
	t150 = t103 * t83;
	t148 = t91 * t111;
	t126 = 0.1e1 + t142;
	t123 = t126 * t93;
	t86 = t102 * t123;
	t147 = qJD(4) * t86;
	t140 = t109 * t111;
	t137 = 0.2e1 * t157;
	t89 = t109 * t153 + 0.1e1;
	t136 = 0.2e1 * (t110 * qJD(4) * t111 * t153 + t109 * t162) / t89 ^ 2;
	t135 = t82 * t158;
	t133 = t107 * t108 * t93;
	t131 = t111 * t159;
	t128 = 0.1e1 + t153;
	t127 = 0.2e1 * t83 * t158;
	t125 = -0.2e1 * t108 * t155;
	t124 = t102 * t133;
	t121 = t128 * t140;
	t87 = 0.1e1 / t89;
	t75 = (-t91 * t124 + t90 * t131) * t103;
	t74 = -t93 * t141 * t146 + (qJD(4) * t123 + t111 * t125) * t103;
	t73 = t86 * t149 + t148 + (-t86 * t148 - t149) * t102;
	t72 = t123 * t144 + 0.2e1 * (t120 * t93 - t126 * t155) * t102;
	t70 = -t87 * qJD(4) * t121 + (t128 * t136 - 0.2e1 * t87 * t162) * t108;
	t68 = (-t138 * t154 + (0.2e1 * t135 + (t105 * t75 + t71) * t156) * t111) * t102 + (t75 * t127 * t111 + (-t75 * t130 + (t75 * t137 + ((-t76 * t124 - t159 * t138 + t155 * t160) * t90 + (-t76 * t131 + (t107 * t125 + (t106 * t109 + t160) * t93 * qJD(4)) * t102) * t91) * t150) * t111 + (-t82 + (-(t97 - t98) * t91 * t133 + t159 * t132) * t83) * t143) * t78) * t103;
	t1 = [t74, 0, t74, t72, 0; t68, 0, t68, (-t146 * t154 + (-0.2e1 * t135 + (-qJD(4) * t73 - t71) * t156) * t103) * t112 + (t73 * t103 * t127 + (-t103 * qJD(4) * t82 + (t103 * t137 + t134) * t73 + (-((t72 - t144) * t90 + (t76 * t86 + qJD(4) + (-t76 - t147) * t102) * t91) * t112 + (-(-t102 * t72 - t144 * t86) * t91 - (t152 * t86 - t147 - t163) * t90) * t111) * t150) * t78) * t111, 0; t70, 0, t70, t136 * t140 * t151 + (-t105 * t121 + (-0.2e1 * t107 * t110 - t108) * t99 * t139) * t87, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end