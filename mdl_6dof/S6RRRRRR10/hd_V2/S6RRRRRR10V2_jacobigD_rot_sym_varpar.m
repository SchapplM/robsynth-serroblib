% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10V2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)), 0, 0, 0, 0; 0, qJD(1) * sin(qJ(1)), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = qJD(1) * cos(qJ(1));
	t27 = qJD(1) * sin(qJ(1));
	t1 = [0, t28, t28, 0, 0, 0; 0, t27, t27, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->4), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
	t115 = qJD(2) + qJD(3);
	t116 = qJ(2) + qJ(3);
	t119 = cos(t116) * t115;
	t117 = sin(qJ(1));
	t111 = qJD(1) * t117;
	t118 = cos(qJ(1));
	t112 = qJD(1) * t118;
	t113 = sin(t116);
	t1 = [0, t112, t112, -t113 * t111 + t118 * t119, 0, 0; 0, t111, t111, t113 * t112 + t117 * t119, 0, 0; 0, 0, 0, t115 * t113, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (30->11), mult. (43->24), div. (0->0), fcn. (43->6), ass. (0->15)
	t157 = qJD(2) + qJD(3);
	t160 = sin(qJ(1));
	t166 = t157 * t160;
	t162 = cos(qJ(1));
	t165 = t157 * t162;
	t153 = qJD(1) * t160;
	t154 = qJD(1) * t162;
	t158 = qJ(2) + qJ(3);
	t156 = cos(t158);
	t164 = qJD(1) * t156 - qJD(4);
	t161 = cos(qJ(4));
	t163 = (qJD(4) * t156 - qJD(1)) * t161;
	t159 = sin(qJ(4));
	t155 = sin(t158);
	t1 = [0, t154, t154, -t155 * t153 + t156 * t165, t162 * t163 + (-t155 * t165 - t164 * t160) * t159, 0; 0, t153, t153, t155 * t154 + t156 * t166, t160 * t163 + (-t155 * t166 + t164 * t162) * t159, 0; 0, 0, 0, t157 * t155, t155 * qJD(4) * t161 + t157 * t156 * t159, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:03
	% EndTime: 2019-10-10 13:38:03
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (77->28), mult. (112->57), div. (0->0), fcn. (114->8), ass. (0->26)
	t232 = qJ(2) + qJ(3);
	t230 = cos(t232);
	t237 = cos(qJ(4));
	t251 = t230 * t237;
	t231 = qJD(2) + qJD(3);
	t250 = t231 * t230;
	t235 = sin(qJ(1));
	t249 = t231 * t235;
	t248 = t231 * t237;
	t233 = sin(qJ(5));
	t247 = t233 * t237;
	t238 = cos(qJ(1));
	t246 = t238 * t231;
	t227 = qJD(1) * t235;
	t228 = qJD(1) * t238;
	t234 = sin(qJ(4));
	t245 = qJD(4) * t234;
	t244 = -qJD(5) + t248;
	t243 = qJD(4) * t230 - qJD(1);
	t242 = qJD(1) * t230 - qJD(4);
	t241 = t243 * t237;
	t240 = t242 * t235;
	t229 = sin(t232);
	t239 = -t229 * t227 + t230 * t246;
	t236 = cos(qJ(5));
	t1 = [0, t228, t228, t239, t238 * t241 + (-t229 * t246 - t240) * t234, -t240 * t247 + (-t244 * t229 - t243 * t234) * t233 * t238 + ((t235 * t234 + t238 * t251) * qJD(5) - t239) * t236; 0, t227, t227, t229 * t228 + t230 * t249, t235 * t241 + (-t229 * t249 + t242 * t238) * t234, (t242 * t247 + (-qJD(1) * t229 - t234 * qJD(5)) * t236) * t238 + ((qJD(1) * t234 - t229 * t248 - t230 * t245) * t233 - t236 * t250 + (t229 * t233 + t236 * t251) * qJD(5)) * t235; 0, 0, 0, t231 * t229, t229 * qJD(4) * t237 + t234 * t250, t244 * t233 * t230 + (-t233 * t245 + (qJD(5) * t237 - t231) * t236) * t229;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end