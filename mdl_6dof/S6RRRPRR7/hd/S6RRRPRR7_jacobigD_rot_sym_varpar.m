% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR7_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t109 = sin(qJ(2));
	t110 = sin(qJ(1));
	t117 = t109 * t110;
	t112 = cos(qJ(1));
	t116 = t109 * t112;
	t111 = cos(qJ(2));
	t115 = t110 * t111;
	t114 = t111 * t112;
	t107 = sin(pkin(6));
	t113 = qJD(1) * t107;
	t108 = cos(pkin(6));
	t1 = [0, t112 * t113, (-t108 * t117 + t114) * qJD(2) + (t108 * t114 - t117) * qJD(1), 0, 0, 0; 0, t110 * t113, (t108 * t116 + t115) * qJD(2) + (t108 * t115 + t116) * qJD(1), 0, 0, 0; 0, 0, t107 * qJD(2) * t109, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t131 = sin(qJ(2));
	t132 = sin(qJ(1));
	t139 = t131 * t132;
	t134 = cos(qJ(1));
	t138 = t131 * t134;
	t133 = cos(qJ(2));
	t137 = t132 * t133;
	t136 = t133 * t134;
	t129 = sin(pkin(6));
	t135 = qJD(1) * t129;
	t130 = cos(pkin(6));
	t128 = t129 * qJD(2) * t131;
	t127 = (t130 * t138 + t137) * qJD(2) + (t130 * t137 + t138) * qJD(1);
	t126 = (-t130 * t139 + t136) * qJD(2) + (t130 * t136 - t139) * qJD(1);
	t1 = [0, t134 * t135, t126, 0, t126, 0; 0, t132 * t135, t127, 0, t127, 0; 0, 0, t128, 0, t128, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:13
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (56->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
	t200 = sin(pkin(6));
	t203 = sin(qJ(1));
	t218 = t200 * t203;
	t205 = cos(qJ(1));
	t217 = t200 * t205;
	t202 = sin(qJ(2));
	t216 = t202 * t203;
	t215 = t202 * t205;
	t204 = cos(qJ(2));
	t214 = t203 * t204;
	t213 = t205 * t204;
	t212 = qJD(1) * t200;
	t198 = qJ(3) + pkin(12) + qJ(5);
	t196 = sin(t198);
	t211 = qJD(2) * t196;
	t210 = qJD(2) * t200;
	t201 = cos(pkin(6));
	t209 = t201 * t213 - t216;
	t208 = t201 * t214 + t215;
	t207 = t201 * t215 + t214;
	t206 = -t201 * t216 + t213;
	t199 = qJD(3) + qJD(5);
	t197 = cos(t198);
	t195 = t202 * t210;
	t194 = qJD(1) * t208 + qJD(2) * t207;
	t193 = qJD(1) * t209 + qJD(2) * t206;
	t1 = [0, t205 * t212, t193, 0, t193, (t196 * t218 + t197 * t206) * t199 - t208 * t211 + (-t196 * t207 - t197 * t217) * qJD(1); 0, t203 * t212, t194, 0, t194, (-t196 * t217 + t197 * t207) * t199 + t209 * t211 + (t196 * t206 - t197 * t218) * qJD(1); 0, 0, t195, 0, t195, t200 * t202 * t199 * t197 + (t199 * t201 + t204 * t210) * t196;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end