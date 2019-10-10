% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(9);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(9);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (28->9), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t167 = sin(qJ(3));
	t172 = qJD(1) * t167;
	t168 = cos(qJ(3));
	t171 = qJD(1) * t168;
	t170 = qJD(3) * t167;
	t169 = qJD(3) * t168;
	t166 = qJ(1) + pkin(9);
	t165 = cos(t166);
	t164 = sin(t166);
	t163 = -t164 * t170 + t165 * t171;
	t162 = -t164 * t169 - t165 * t172;
	t161 = -t164 * t171 - t165 * t170;
	t160 = t164 * t172 - t165 * t169;
	t1 = [-t163, 0, t160, 0, 0, 0; t161, 0, t162, 0, 0, 0; 0, 0, -t170, 0, 0, 0; -qJD(1) * t164, 0, 0, 0, 0, 0; qJD(1) * t165, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t162, 0, t161, 0, 0, 0; -t160, 0, t163, 0, 0, 0; 0, 0, t169, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t50 = sin(qJ(3));
	t55 = qJD(1) * t50;
	t51 = cos(qJ(3));
	t54 = qJD(1) * t51;
	t53 = qJD(3) * t50;
	t52 = qJD(3) * t51;
	t49 = qJ(1) + pkin(9);
	t48 = cos(t49);
	t47 = sin(t49);
	t46 = -t47 * t53 + t48 * t54;
	t45 = t47 * t52 + t48 * t55;
	t44 = t47 * t54 + t48 * t53;
	t43 = -t47 * t55 + t48 * t52;
	t1 = [-t45, 0, -t44, 0, 0, 0; t43, 0, t46, 0, 0, 0; 0, 0, t52, 0, 0, 0; t46, 0, t43, 0, 0, 0; t44, 0, t45, 0, 0, 0; 0, 0, t53, 0, 0, 0; qJD(1) * t47, 0, 0, 0, 0, 0; -qJD(1) * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:45
	% EndTime: 2019-10-10 00:18:45
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (109->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t249 = cos(qJ(6));
	t248 = sin(qJ(3));
	t266 = qJD(1) * t248;
	t255 = qJD(6) + t266;
	t268 = t249 * t255;
	t262 = qJD(6) * t248;
	t256 = qJD(1) + t262;
	t247 = sin(qJ(6));
	t250 = cos(qJ(3));
	t263 = qJD(3) * t250;
	t258 = t247 * t263;
	t267 = t256 * t249 + t258;
	t265 = qJD(1) * t250;
	t264 = qJD(3) * t248;
	t261 = qJD(6) * t250;
	t260 = t247 * t265;
	t259 = t249 * t265;
	t257 = t249 * t263;
	t254 = t255 * t247;
	t253 = -t247 * t261 - t249 * t264;
	t252 = t247 * t264 - t249 * t261;
	t251 = t256 * t247 - t257;
	t246 = qJ(1) + pkin(9);
	t245 = cos(t246);
	t244 = sin(t246);
	t243 = t251 * t244 - t245 * t268;
	t242 = t267 * t244 + t245 * t254;
	t241 = t244 * t268 + t251 * t245;
	t240 = t244 * t254 - t267 * t245;
	t1 = [t243, 0, -t244 * t259 + t253 * t245, 0, 0, t240; -t241, 0, t253 * t244 + t245 * t259, 0, 0, -t242; 0, 0, -t247 * t262 + t257, 0, 0, -t252; t242, 0, t244 * t260 + t252 * t245, 0, 0, t241; t240, 0, t252 * t244 - t245 * t260, 0, 0, t243; 0, 0, -t249 * t262 - t258, 0, 0, t253; t244 * t264 - t245 * t265, 0, t244 * t266 - t245 * t263, 0, 0, 0; -t244 * t265 - t245 * t264, 0, -t244 * t263 - t245 * t266, 0, 0, 0; 0, 0, -t264, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end