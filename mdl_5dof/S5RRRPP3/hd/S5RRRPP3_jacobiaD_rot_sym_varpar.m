% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP3
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
%   Wie in S5RRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:01
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:06
	% EndTime: 2019-12-29 19:40:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:06
	% EndTime: 2019-12-29 19:40:08
	% DurationCPUTime: 1.36s
	% Computational Cost: add. (2226->70), mult. (2659->154), div. (622->14), fcn. (3165->7), ass. (0->75)
	t113 = qJ(1) + qJ(2);
	t104 = sin(t113);
	t105 = cos(t113);
	t152 = t105 / t104;
	t141 = qJD(3) * t104;
	t115 = cos(qJ(3));
	t110 = 0.1e1 / t115;
	t140 = qJD(3) * t115;
	t114 = sin(qJ(3));
	t108 = t114 ^ 2;
	t111 = 0.1e1 / t115 ^ 2;
	t143 = t108 * t111;
	t106 = qJD(1) + qJD(2);
	t144 = t106 * t114;
	t98 = t104 ^ 2;
	t97 = t98 * t143 + 0.1e1;
	t95 = 0.1e1 / t97;
	t78 = (-(-t104 * t140 - t105 * t144) * t110 + t141 * t143) * t95;
	t163 = t78 - t141;
	t103 = t105 ^ 2;
	t162 = t106 * (0.1e1 / t98 * t103 + 0.1e1) * t152;
	t146 = t104 * t114;
	t94 = atan2(-t146, -t115);
	t89 = sin(t94);
	t135 = t89 * t146;
	t90 = cos(t94);
	t87 = -t115 * t90 - t135;
	t84 = 0.1e1 / t87;
	t85 = 0.1e1 / t87 ^ 2;
	t160 = 0.2e1 * t114;
	t159 = t95 - 0.1e1;
	t145 = t105 * t106;
	t126 = t104 * t108 * t145;
	t132 = t114 * t140;
	t151 = t115 * t89;
	t154 = t104 * t78;
	t73 = t163 * t151 + (-t89 * t145 + (qJD(3) - t154) * t90) * t114;
	t157 = t73 * t84 * t85;
	t81 = t103 * t108 * t85 + 0.1e1;
	t158 = (-t85 * t126 + (-t108 * t157 + t85 * t132) * t103) / t81 ^ 2;
	t79 = 0.1e1 / t81;
	t156 = t79 * t85;
	t107 = t114 * t108;
	t109 = t115 ^ 2;
	t142 = t110 * t114;
	t122 = qJD(3) * (t107 * t110 / t109 + t142);
	t155 = (t111 * t126 + t98 * t122) / t97 ^ 2;
	t153 = t105 * t85;
	t150 = t90 * t114;
	t129 = 0.1e1 + t143;
	t124 = t129 * t95;
	t88 = t104 * t124;
	t149 = qJD(3) * t88;
	t148 = 0.1e1 / t104 ^ 2 * t103;
	t147 = t104 * t106;
	t139 = 0.2e1 * t157;
	t93 = t109 * t148 + 0.1e1;
	t138 = 0.2e1 * (-t109 * t162 - t132 * t148) / t93 ^ 2;
	t137 = t84 * t158;
	t136 = t108 * t110 * t95;
	t134 = t114 * t159;
	t133 = t79 * t140;
	t131 = 0.2e1 * t85 * t158;
	t130 = 0.1e1 + t148;
	t128 = -0.2e1 * t110 * t155;
	t127 = t104 * t136;
	t91 = 0.1e1 / t93;
	t125 = t130 * t91;
	t77 = (-t90 * t127 + t89 * t134) * t105;
	t76 = -t95 * t142 * t147 + (qJD(3) * t124 + t114 * t128) * t105;
	t75 = t88 * t151 + t150 + (-t88 * t150 - t151) * t104;
	t74 = t124 * t145 + 0.2e1 * (t122 * t95 - t129 * t155) * t104;
	t72 = t114 * qJD(3) * t125 + (t130 * t138 + 0.2e1 * t91 * t162) * t115;
	t70 = (-t84 * t133 + (0.2e1 * t137 + (t106 * t77 + t73) * t156) * t114) * t104 + (-t77 * t85 * t133 + (t77 * t131 + (t77 * t139 + ((-t78 * t127 - t159 * t140 + t155 * t160) * t89 + (-t78 * t134 + (t108 * t128 + (t107 * t111 + t160) * t95 * qJD(3)) * t104) * t90) * t153) * t79) * t114 + (-t84 + (-(-t103 + t98) * t90 * t136 + t159 * t135) * t85) * t79 * t144) * t105;
	t1 = [t76, t76, t74, 0, 0; t70, t70, (-t84 * t79 * t147 + (-0.2e1 * t137 + (-qJD(3) * t75 - t73) * t156) * t105) * t115 + (t75 * t105 * t131 + (-t105 * qJD(3) * t84 + (t105 * t139 + t85 * t147) * t75 + (-((t74 - t145) * t89 + (t78 * t88 + qJD(3) + (-t78 - t149) * t104) * t90) * t115 + (-(-t104 * t74 - t145 * t88) * t90 - (t154 * t88 - t149 - t163) * t89) * t114) * t153) * t79) * t114, 0, 0; t72, t72, -t91 * t140 * t152 + (t106 * t125 + t138 * t152) * t114, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:40:01
	% EndTime: 2019-12-29 19:40:02
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (1944->69), mult. (2659->154), div. (622->14), fcn. (3165->7), ass. (0->76)
	t115 = sin(qJ(3));
	t108 = t115 ^ 2;
	t114 = qJ(1) + qJ(2);
	t106 = cos(t114);
	t104 = t106 ^ 2;
	t107 = qJD(1) + qJD(2);
	t105 = sin(t114);
	t153 = 0.1e1 / t105 * t106;
	t99 = t105 ^ 2;
	t124 = t107 * (-0.1e1 / t99 * t104 - 0.1e1) * t153;
	t116 = cos(qJ(3));
	t142 = qJD(3) * t116;
	t134 = t115 * t142;
	t152 = 0.1e1 / t105 ^ 2 * t104;
	t94 = t108 * t152 + 0.1e1;
	t165 = (t108 * t124 + t134 * t152) / t94 ^ 2;
	t110 = 0.1e1 / t115 ^ 2;
	t109 = 0.1e1 / t115;
	t113 = t116 ^ 2;
	t112 = t116 * t113;
	t146 = t109 * t116;
	t123 = qJD(3) * (-t109 / t108 * t112 - t146);
	t149 = t106 * t107;
	t127 = t105 * t113 * t149;
	t145 = t110 * t113;
	t98 = t99 * t145 + 0.1e1;
	t164 = (t110 * t127 + t99 * t123) / t98 ^ 2;
	t143 = qJD(3) * t115;
	t144 = qJD(3) * t105;
	t148 = t107 * t116;
	t96 = 0.1e1 / t98;
	t79 = ((t105 * t143 - t106 * t148) * t109 + t144 * t145) * t96;
	t130 = 0.1e1 + t145;
	t125 = t130 * t96;
	t89 = t105 * t125;
	t163 = qJD(3) * t89 + t79;
	t150 = t105 * t116;
	t95 = atan2(-t150, t115);
	t92 = sin(t95);
	t138 = t92 * t150;
	t93 = cos(t95);
	t88 = t93 * t115 - t138;
	t85 = 0.1e1 / t88;
	t86 = 0.1e1 / t88 ^ 2;
	t161 = t96 - 0.1e1;
	t156 = t115 * t92;
	t158 = t105 * t79;
	t74 = (-t79 + t144) * t156 + (-t92 * t149 + (qJD(3) - t158) * t93) * t116;
	t160 = t74 * t85 * t86;
	t82 = t104 * t113 * t86 + 0.1e1;
	t80 = 0.1e1 / t82;
	t159 = t80 * t86;
	t157 = t106 * t86;
	t155 = t116 * t93;
	t151 = t105 * t107;
	t147 = t109 * t113;
	t141 = 0.2e1 * (-t86 * t127 + (-t113 * t160 - t86 * t134) * t104) / t82 ^ 2;
	t140 = 0.2e1 * t160;
	t139 = 0.2e1 * t164;
	t137 = t96 * t147;
	t136 = t161 * t116;
	t135 = t80 * t143;
	t133 = t85 * t141;
	t132 = t86 * t141;
	t131 = 0.1e1 + t152;
	t129 = t116 * t139;
	t128 = t105 * t137;
	t90 = 0.1e1 / t94;
	t126 = t131 * t90;
	t78 = (t93 * t128 + t92 * t136) * t106;
	t77 = t96 * t146 * t151 + (qJD(3) * t125 + t109 * t129) * t106;
	t76 = -t89 * t156 + t155 + (-t89 * t155 + t156) * t105;
	t75 = t125 * t149 + 0.2e1 * (t123 * t96 - t130 * t164) * t105;
	t73 = t126 * t142 + 0.2e1 * (t124 * t90 - t131 * t165) * t115;
	t71 = (t85 * t135 + (t133 + (t107 * t78 + t74) * t159) * t116) * t105 + (t78 * t86 * t135 + (t78 * t132 + (t78 * t140 + ((t79 * t128 + t161 * t143 + t129) * t92 + (-t79 * t136 + (t139 * t147 + (t110 * t112 + 0.2e1 * t116) * t96 * qJD(3)) * t105) * t93) * t157) * t80) * t116 + (-t85 + (-(t104 - t99) * t93 * t137 + t161 * t138) * t86) * t80 * t148) * t106;
	t1 = [t77, t77, t75, 0, 0; t71, t71, (t85 * t80 * t151 + (t133 + (qJD(3) * t76 + t74) * t159) * t106) * t115 + (t76 * t106 * t132 + (-t106 * qJD(3) * t85 + (t106 * t140 + t86 * t151) * t76 + (-((-t75 + t149) * t92 + (t163 * t105 - t79 * t89 - qJD(3)) * t93) * t115 + (-(-t105 * t75 - t149 * t89) * t93 - (t158 * t89 + t144 - t163) * t92) * t116) * t157) * t80) * t116, 0, 0; t73, t73, t90 * t143 * t153 + (t107 * t126 + 0.2e1 * t153 * t165) * t116, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end