% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRPR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t22 = qJD(1) * sin(qJ(1));
	t21 = qJD(1) * cos(qJ(1));
	t18 = cos(pkin(8));
	t17 = sin(pkin(8));
	t1 = [-t18 * t21, 0, 0, 0, 0; -t18 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t17 * t21, 0, 0, 0, 0; t17 * t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0; t21, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:40
	% EndTime: 2019-12-31 18:33:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t48 = qJD(1) * t43;
	t44 = cos(qJ(1));
	t47 = qJD(1) * t44;
	t46 = qJD(3) * t43;
	t45 = qJD(3) * t44;
	t42 = pkin(8) + qJ(3);
	t41 = cos(t42);
	t40 = sin(t42);
	t39 = t40 * t46 - t41 * t47;
	t38 = t40 * t47 + t41 * t46;
	t37 = t40 * t45 + t41 * t48;
	t36 = t40 * t48 - t41 * t45;
	t1 = [t39, 0, t36, 0, 0; -t37, 0, -t38, 0, 0; 0, 0, -qJD(3) * t40, 0, 0; t38, 0, t37, 0, 0; t36, 0, t39, 0, 0; 0, 0, -qJD(3) * t41, 0, 0; -t48, 0, 0, 0, 0; t47, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:41
	% EndTime: 2019-12-31 18:33:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t170 = sin(qJ(1));
	t175 = qJD(1) * t170;
	t171 = cos(qJ(1));
	t174 = qJD(1) * t171;
	t173 = qJD(3) * t170;
	t172 = qJD(3) * t171;
	t169 = pkin(8) + qJ(3);
	t168 = cos(t169);
	t167 = sin(t169);
	t166 = -t167 * t173 + t168 * t174;
	t165 = t167 * t174 + t168 * t173;
	t164 = t167 * t172 + t168 * t175;
	t163 = -t167 * t175 + t168 * t172;
	t1 = [-t175, 0, 0, 0, 0; t174, 0, 0, 0, 0; 0, 0, 0, 0, 0; t166, 0, t163, 0, 0; t164, 0, t165, 0, 0; 0, 0, qJD(3) * t167, 0, 0; -t165, 0, -t164, 0, 0; t163, 0, t166, 0, 0; 0, 0, qJD(3) * t168, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:33:41
	% EndTime: 2019-12-31 18:33:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->29), mult. (173->56), div. (0->0), fcn. (173->6), ass. (0->32)
	t246 = cos(qJ(5));
	t243 = pkin(8) + qJ(3);
	t241 = sin(t243);
	t251 = qJD(5) * t241 + qJD(1);
	t267 = t246 * t251;
	t244 = sin(qJ(5));
	t266 = t251 * t244;
	t245 = sin(qJ(1));
	t265 = qJD(1) * t245;
	t247 = cos(qJ(1));
	t264 = qJD(1) * t247;
	t263 = qJD(3) * t244;
	t262 = qJD(3) * t245;
	t261 = qJD(3) * t246;
	t260 = qJD(3) * t247;
	t259 = qJD(5) * t244;
	t258 = qJD(5) * t246;
	t257 = qJD(5) * t247;
	t242 = cos(t243);
	t256 = t242 * t261;
	t255 = t241 * t262;
	t254 = t242 * t262;
	t253 = t241 * t260;
	t252 = t242 * t260;
	t250 = -qJD(1) * t241 - qJD(5);
	t249 = t250 * t247;
	t248 = t250 * t245 + t252;
	t240 = t248 * t244 + t247 * t267;
	t239 = t248 * t246 - t247 * t266;
	t238 = -t245 * t267 + (t249 - t254) * t244;
	t237 = t246 * t249 + (-t256 + t266) * t245;
	t1 = [t238, 0, -t244 * t253 + (-t244 * t265 + t246 * t257) * t242, 0, t239; t240, 0, -t244 * t255 + (t244 * t264 + t245 * t258) * t242, 0, -t237; 0, 0, t241 * t258 + t242 * t263, 0, t241 * t261 + t242 * t259; t237, 0, -t246 * t253 + (-t244 * t257 - t246 * t265) * t242, 0, -t240; t239, 0, -t246 * t255 + (-t245 * t259 + t246 * t264) * t242, 0, t238; 0, 0, -t241 * t259 + t256, 0, -t241 * t263 + t242 * t258; -t242 * t264 + t255, 0, t241 * t265 - t252, 0, 0; -t242 * t265 - t253, 0, -t241 * t264 - t254, 0, 0; 0, 0, -qJD(3) * t241, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end