% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(11));
	t79 = sin(pkin(11));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t125 = sin(pkin(6));
	t128 = sin(qJ(3));
	t138 = t125 * t128;
	t129 = sin(qJ(2));
	t137 = t125 * t129;
	t127 = cos(pkin(6));
	t136 = t127 * t129;
	t131 = cos(qJ(2));
	t135 = t127 * t131;
	t134 = qJD(2) * t128;
	t124 = sin(pkin(11));
	t126 = cos(pkin(11));
	t133 = t124 * t131 + t126 * t136;
	t132 = -t124 * t136 + t126 * t131;
	t130 = cos(qJ(3));
	t1 = [0, 0, t132 * qJD(2), (t124 * t138 + t132 * t130) * qJD(3) + (-t124 * t135 - t126 * t129) * t134, 0, 0; 0, 0, t133 * qJD(2), (-t126 * t138 + t133 * t130) * qJD(3) + (-t124 * t129 + t126 * t135) * t134, 0, 0; 0, 0, qJD(2) * t137, t125 * t131 * t134 + (t127 * t128 + t130 * t137) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t153 = sin(pkin(6));
	t156 = sin(qJ(3));
	t166 = t153 * t156;
	t157 = sin(qJ(2));
	t165 = t153 * t157;
	t155 = cos(pkin(6));
	t164 = t155 * t157;
	t159 = cos(qJ(2));
	t163 = t155 * t159;
	t162 = qJD(2) * t156;
	t152 = sin(pkin(11));
	t154 = cos(pkin(11));
	t161 = t152 * t159 + t154 * t164;
	t160 = -t152 * t164 + t154 * t159;
	t158 = cos(qJ(3));
	t1 = [0, 0, t160 * qJD(2), (t152 * t166 + t160 * t158) * qJD(3) + (-t152 * t163 - t154 * t157) * t162, 0, 0; 0, 0, t161 * qJD(2), (-t154 * t166 + t161 * t158) * qJD(3) + (-t152 * t157 + t154 * t163) * t162, 0, 0; 0, 0, qJD(2) * t165, t153 * t159 * t162 + (t155 * t156 + t158 * t165) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:08
	% EndTime: 2019-10-09 22:56:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->13), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
	t183 = sin(pkin(6));
	t186 = sin(qJ(3));
	t196 = t183 * t186;
	t187 = sin(qJ(2));
	t195 = t183 * t187;
	t185 = cos(pkin(6));
	t194 = t185 * t187;
	t189 = cos(qJ(2));
	t193 = t185 * t189;
	t192 = qJD(2) * t186;
	t182 = sin(pkin(11));
	t184 = cos(pkin(11));
	t191 = t182 * t189 + t184 * t194;
	t190 = -t182 * t194 + t184 * t189;
	t188 = cos(qJ(3));
	t181 = t183 * t189 * t192 + (t185 * t186 + t188 * t195) * qJD(3);
	t180 = (t182 * t196 + t190 * t188) * qJD(3) + (-t182 * t193 - t184 * t187) * t192;
	t179 = (-t184 * t196 + t191 * t188) * qJD(3) + (-t182 * t187 + t184 * t193) * t192;
	t1 = [0, 0, t190 * qJD(2), t180, 0, -t180; 0, 0, t191 * qJD(2), t179, 0, -t179; 0, 0, qJD(2) * t195, t181, 0, -t181;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end