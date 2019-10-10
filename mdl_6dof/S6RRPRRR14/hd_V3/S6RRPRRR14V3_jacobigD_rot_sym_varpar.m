% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14V3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobigD_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (11->8), div. (0->0), fcn. (11->4), ass. (0->7)
	t72 = sin(qJ(1));
	t77 = qJD(1) * t72;
	t74 = cos(qJ(1));
	t76 = qJD(1) * t74;
	t75 = qJD(2) * cos(qJ(2));
	t71 = sin(qJ(2));
	t1 = [0, t76, 0, -t71 * t77 + t74 * t75, 0, 0; 0, t77, 0, t71 * t76 + t72 * t75, 0, 0; 0, 0, 0, qJD(2) * t71, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (41->23), div. (0->0), fcn. (41->6), ass. (0->14)
	t116 = sin(qJ(1));
	t126 = qJD(1) * t116;
	t119 = cos(qJ(1));
	t125 = qJD(1) * t119;
	t115 = sin(qJ(2));
	t124 = qJD(2) * t115;
	t118 = cos(qJ(2));
	t123 = qJD(2) * t118;
	t122 = qJD(2) * t119;
	t121 = qJD(1) * t118 - qJD(4);
	t117 = cos(qJ(4));
	t120 = (qJD(4) * t118 - qJD(1)) * t117;
	t114 = sin(qJ(4));
	t1 = [0, t125, 0, -t115 * t126 + t118 * t122, t119 * t120 + (-t115 * t122 - t121 * t116) * t114, 0; 0, t126, 0, t115 * t125 + t116 * t123, t116 * t120 + (-t116 * t124 + t121 * t119) * t114, 0; 0, 0, 0, t124, t115 * qJD(4) * t117 + t114 * t123, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:24
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (33->25), mult. (110->53), div. (0->0), fcn. (112->8), ass. (0->25)
	t181 = sin(qJ(4));
	t186 = cos(qJ(2));
	t204 = -qJD(4) * t186 + qJD(1);
	t205 = t204 * t181;
	t180 = sin(qJ(5));
	t185 = cos(qJ(4));
	t203 = t180 * t185;
	t202 = t185 * t186;
	t183 = sin(qJ(1));
	t201 = qJD(1) * t183;
	t187 = cos(qJ(1));
	t200 = qJD(1) * t187;
	t182 = sin(qJ(2));
	t199 = qJD(2) * t182;
	t198 = qJD(2) * t185;
	t197 = qJD(2) * t186;
	t196 = qJD(4) * t182;
	t194 = t187 * qJD(2);
	t192 = qJD(1) * t186 - qJD(4);
	t191 = -qJD(5) + t198;
	t190 = t204 * t185;
	t189 = t192 * t183;
	t188 = -t182 * t201 + t186 * t194;
	t184 = cos(qJ(5));
	t1 = [0, t200, 0, t188, -t187 * t190 + (-t182 * t194 - t189) * t181, -t189 * t203 + (-t191 * t182 + t205) * t180 * t187 + ((t183 * t181 + t187 * t202) * qJD(5) - t188) * t184; 0, t201, 0, t182 * t200 + t183 * t197, -t183 * t190 + (-t183 * t199 + t192 * t187) * t181, (t192 * t203 + (-qJD(1) * t182 - t181 * qJD(5)) * t184) * t187 + ((-t182 * t198 + t205) * t180 - t184 * t197 + (t182 * t180 + t184 * t202) * qJD(5)) * t183; 0, 0, 0, t199, t181 * t197 + t185 * t196, (qJD(5) * t185 - qJD(2)) * t184 * t182 + (-t181 * t196 + t191 * t186) * t180;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end