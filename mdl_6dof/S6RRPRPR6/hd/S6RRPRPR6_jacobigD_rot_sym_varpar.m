% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t133 = sin(pkin(11));
	t135 = cos(pkin(11));
	t137 = sin(qJ(2));
	t139 = cos(qJ(2));
	t142 = t139 * t133 + t137 * t135;
	t144 = qJD(2) * t142;
	t134 = sin(pkin(6));
	t143 = qJD(1) * t134;
	t141 = t133 * t137 - t135 * t139;
	t140 = cos(qJ(1));
	t138 = sin(qJ(1));
	t136 = cos(pkin(6));
	t131 = t141 * qJD(2);
	t130 = t141 * t136;
	t129 = t136 * t144;
	t1 = [0, t140 * t143, 0, -t138 * t129 - t140 * t131 + (-t130 * t140 - t138 * t142) * qJD(1), 0, 0; 0, t138 * t143, 0, t140 * t129 - t138 * t131 + (-t130 * t138 + t140 * t142) * qJD(1), 0, 0; 0, 0, 0, t134 * t144, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t145 = sin(pkin(11));
	t147 = cos(pkin(11));
	t149 = sin(qJ(2));
	t151 = cos(qJ(2));
	t154 = t151 * t145 + t149 * t147;
	t156 = qJD(2) * t154;
	t146 = sin(pkin(6));
	t155 = qJD(1) * t146;
	t153 = t145 * t149 - t147 * t151;
	t152 = cos(qJ(1));
	t150 = sin(qJ(1));
	t148 = cos(pkin(6));
	t143 = t153 * qJD(2);
	t142 = t153 * t148;
	t141 = t148 * t156;
	t1 = [0, t152 * t155, 0, -t150 * t141 - t152 * t143 + (-t142 * t152 - t150 * t154) * qJD(1), 0, 0; 0, t150 * t155, 0, t152 * t141 - t150 * t143 + (-t142 * t150 + t152 * t154) * qJD(1), 0, 0; 0, 0, 0, t146 * t156, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
	t197 = sin(pkin(11));
	t199 = cos(pkin(11));
	t202 = sin(qJ(2));
	t205 = cos(qJ(2));
	t208 = t202 * t197 - t205 * t199;
	t215 = t208 * qJD(2);
	t209 = t205 * t197 + t202 * t199;
	t194 = t209 * qJD(2);
	t198 = sin(pkin(6));
	t203 = sin(qJ(1));
	t214 = t198 * t203;
	t206 = cos(qJ(1));
	t213 = t198 * t206;
	t212 = qJD(1) * t198;
	t200 = cos(pkin(6));
	t192 = t209 * t200;
	t211 = t206 * t192 - t203 * t208;
	t210 = -t203 * t192 - t206 * t208;
	t204 = cos(qJ(4));
	t201 = sin(qJ(4));
	t191 = t208 * t200;
	t190 = t200 * t215;
	t189 = t200 * t194;
	t1 = [0, t206 * t212, 0, -t203 * t189 - t206 * t215 + (-t191 * t206 - t203 * t209) * qJD(1), 0, (t203 * t190 - t206 * t194) * t204 + (-t210 * t201 + t204 * t214) * qJD(4) + (t201 * t213 - t211 * t204) * qJD(1); 0, t203 * t212, 0, t206 * t189 - t203 * t215 + (-t191 * t203 + t206 * t209) * qJD(1), 0, (-t206 * t190 - t203 * t194) * t204 + (-t211 * t201 - t204 * t213) * qJD(4) + (t201 * t214 + t210 * t204) * qJD(1); 0, 0, 0, t198 * t194, 0, t200 * qJD(4) * t204 + (-t209 * qJD(4) * t201 - t204 * t215) * t198;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end