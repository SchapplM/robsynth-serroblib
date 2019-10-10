% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:24
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->10), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t157 = sin(qJ(1));
	t164 = qJD(1) * t157;
	t159 = cos(qJ(1));
	t163 = qJD(1) * t159;
	t156 = sin(qJ(3));
	t162 = qJD(3) * t156;
	t158 = cos(qJ(3));
	t161 = qJD(3) * t158;
	t160 = qJD(3) * t159;
	t155 = -t157 * t162 + t158 * t163;
	t154 = t156 * t163 + t157 * t161;
	t153 = t156 * t160 + t158 * t164;
	t152 = t156 * t164 - t158 * t160;
	t1 = [-t152, 0, t155, 0, 0, 0; t154, 0, t153, 0, 0, 0; 0, 0, -t161, 0, 0, 0; -t163, 0, 0, 0, 0, 0; -t164, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t153, 0, t154, 0, 0, 0; -t155, 0, t152, 0, 0, 0; 0, 0, -t162, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t43 = sin(qJ(1));
	t50 = qJD(1) * t43;
	t45 = cos(qJ(1));
	t49 = qJD(1) * t45;
	t42 = sin(qJ(3));
	t48 = qJD(3) * t42;
	t44 = cos(qJ(3));
	t47 = qJD(3) * t44;
	t46 = qJD(3) * t45;
	t41 = t43 * t48 - t44 * t49;
	t40 = t42 * t49 + t43 * t47;
	t39 = t42 * t46 + t44 * t50;
	t38 = t42 * t50 - t44 * t46;
	t1 = [t39, 0, t40, 0, 0, 0; t41, 0, t38, 0, 0, 0; 0, 0, -t48, 0, 0, 0; t38, 0, t41, 0, 0, 0; -t40, 0, -t39, 0, 0, 0; 0, 0, t47, 0, 0, 0; t49, 0, 0, 0, 0, 0; t50, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:24
	% EndTime: 2019-10-10 00:27:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->26), mult. (173->49), div. (0->0), fcn. (173->6), ass. (0->34)
	t239 = cos(qJ(6));
	t240 = cos(qJ(3));
	t254 = qJD(6) * t240;
	t247 = qJD(1) + t254;
	t264 = t239 * t247;
	t236 = sin(qJ(6));
	t263 = t247 * t236;
	t238 = sin(qJ(1));
	t262 = qJD(1) * t238;
	t261 = qJD(1) * t240;
	t241 = cos(qJ(1));
	t260 = qJD(1) * t241;
	t237 = sin(qJ(3));
	t259 = qJD(3) * t237;
	t258 = qJD(3) * t240;
	t257 = qJD(3) * t241;
	t256 = qJD(6) * t236;
	t255 = qJD(6) * t237;
	t253 = t239 * t259;
	t252 = t239 * t255;
	t251 = t238 * t259;
	t250 = t238 * t258;
	t249 = t237 * t257;
	t248 = t240 * t257;
	t246 = qJD(6) + t261;
	t245 = t246 * t241;
	t244 = t237 * t260 + t250;
	t243 = -t237 * t262 + t248;
	t242 = t246 * t238 + t249;
	t235 = t239 * t245 + (-t253 - t263) * t238;
	t234 = t238 * t264 + (t245 - t251) * t236;
	t233 = t242 * t239 + t241 * t263;
	t232 = t242 * t236 - t241 * t264;
	t1 = [t233, 0, t239 * t250 + (-t238 * t256 + t239 * t260) * t237, 0, 0, t234; -t235, 0, -t239 * t248 + (t239 * t262 + t241 * t256) * t237, 0, 0, t232; 0, 0, -t236 * t254 - t253, 0, 0, -t236 * t258 - t252; -t232, 0, -t244 * t236 - t238 * t252, 0, 0, t235; t234, 0, t243 * t236 + t241 * t252, 0, 0, t233; 0, 0, t236 * t259 - t239 * t254, 0, 0, t236 * t255 - t239 * t258; t243, 0, t240 * t260 - t251, 0, 0, 0; t244, 0, t238 * t261 + t249, 0, 0, 0; 0, 0, -t258, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end