% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRRR15_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
JgD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(5));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (8->3), div. (0->0), fcn. (8->3), ass. (0->4)
	t101 = qJD(1) * sin(pkin(5));
	t99 = cos(qJ(1)) * t101;
	t98 = sin(qJ(1)) * t101;
	t1 = [0, t99, t99, 0, 0; 0, t98, t98, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (12->3), div. (0->0), fcn. (12->3), ass. (0->4)
	t109 = qJD(1) * sin(pkin(5));
	t107 = cos(qJ(1)) * t109;
	t106 = sin(qJ(1)) * t109;
	t1 = [0, t107, t107, t107, 0; 0, t106, t106, t106, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (36->10), mult. (49->24), div. (0->0), fcn. (49->8), ass. (0->16)
	t221 = cos(pkin(5));
	t222 = sin(qJ(1));
	t226 = t221 * t222;
	t223 = cos(qJ(1));
	t225 = t221 * t223;
	t219 = sin(pkin(5));
	t224 = qJD(1) * t219;
	t212 = t222 * t224;
	t213 = t223 * t224;
	t220 = cos(pkin(6));
	t218 = sin(pkin(6));
	t217 = qJ(2) + qJ(3) + qJ(4);
	t216 = qJD(2) + qJD(3) + qJD(4);
	t215 = cos(t217);
	t214 = sin(t217);
	t1 = [0, t213, t213, t213, t220 * t213 + ((-t214 * t226 + t215 * t223) * t216 + (-t214 * t222 + t215 * t225) * qJD(1)) * t218; 0, t212, t212, t212, t220 * t212 + ((t214 * t225 + t215 * t222) * t216 + (t214 * t223 + t215 * t226) * qJD(1)) * t218; 0, 0, 0, 0, t219 * t216 * t214 * t218;];
	JgD_rot = t1;
end