% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t95 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t95, 0, 0, 0, 0; 0, sin(qJ(1)) * t95, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t139 = sin(pkin(6));
	t148 = qJD(1) * t139;
	t138 = sin(pkin(11));
	t140 = cos(pkin(11));
	t142 = sin(qJ(2));
	t144 = cos(qJ(2));
	t147 = t138 * t144 + t140 * t142;
	t137 = -t142 * t138 + t144 * t140;
	t146 = qJD(2) * t137;
	t145 = cos(qJ(1));
	t143 = sin(qJ(1));
	t141 = cos(pkin(6));
	t136 = t147 * qJD(2);
	t135 = t147 * t141;
	t134 = t141 * t146;
	t1 = [0, t145 * t148, 0, 0, -t143 * t134 - t145 * t136 + (-t135 * t145 - t137 * t143) * qJD(1), 0; 0, t143 * t148, 0, 0, t145 * t134 - t143 * t136 + (-t135 * t143 + t137 * t145) * qJD(1), 0; 0, 0, 0, 0, t139 * t146, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (45->23), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
	t203 = sin(pkin(11));
	t205 = cos(pkin(11));
	t208 = sin(qJ(2));
	t211 = cos(qJ(2));
	t214 = t208 * t203 - t211 * t205;
	t221 = t214 * qJD(2);
	t215 = t211 * t203 + t208 * t205;
	t200 = t215 * qJD(2);
	t204 = sin(pkin(6));
	t209 = sin(qJ(1));
	t220 = t204 * t209;
	t212 = cos(qJ(1));
	t219 = t204 * t212;
	t218 = qJD(1) * t204;
	t206 = cos(pkin(6));
	t197 = t214 * t206;
	t217 = -t212 * t197 - t209 * t215;
	t216 = -t209 * t197 + t212 * t215;
	t210 = cos(qJ(5));
	t207 = sin(qJ(5));
	t198 = t215 * t206;
	t196 = t206 * t221;
	t195 = t206 * t200;
	t1 = [0, t212 * t218, 0, 0, t209 * t196 - t212 * t200 + (-t198 * t212 + t209 * t214) * qJD(1), -(-t209 * t195 - t212 * t221) * t210 + (t216 * t207 + t210 * t220) * qJD(5) + (t207 * t219 - t217 * t210) * qJD(1); 0, t209 * t218, 0, 0, -t212 * t196 - t209 * t200 + (-t198 * t209 - t212 * t214) * qJD(1), -(t212 * t195 - t209 * t221) * t210 + (-t217 * t207 - t210 * t219) * qJD(5) + (t207 * t220 - t216 * t210) * qJD(1); 0, 0, 0, 0, -t204 * t221, t206 * qJD(5) * t210 + (t214 * qJD(5) * t207 - t210 * t200) * t204;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end