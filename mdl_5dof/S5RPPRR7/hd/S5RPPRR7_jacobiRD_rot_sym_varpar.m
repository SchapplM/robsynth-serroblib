% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:07
	% EndTime: 2019-12-29 16:16:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:07
	% EndTime: 2019-12-29 16:16:07
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:07
	% EndTime: 2019-12-29 16:16:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:07
	% EndTime: 2019-12-29 16:16:07
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(8);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:07
	% EndTime: 2019-12-29 16:16:07
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t42 = sin(qJ(4));
	t47 = qJD(1) * t42;
	t43 = cos(qJ(4));
	t46 = qJD(1) * t43;
	t45 = qJD(4) * t42;
	t44 = qJD(4) * t43;
	t41 = qJ(1) + pkin(8);
	t40 = cos(t41);
	t39 = sin(t41);
	t38 = -t39 * t45 + t40 * t46;
	t37 = t39 * t44 + t40 * t47;
	t36 = t39 * t46 + t40 * t45;
	t35 = -t39 * t47 + t40 * t44;
	t1 = [t35, 0, 0, t38, 0; t37, 0, 0, t36, 0; 0, 0, 0, -t44, 0; -t36, 0, 0, -t37, 0; t38, 0, 0, t35, 0; 0, 0, 0, t45, 0; -qJD(1) * t40, 0, 0, 0, 0; -qJD(1) * t39, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:16:09
	% EndTime: 2019-12-29 16:16:09
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (109->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t239 = sin(qJ(5));
	t240 = sin(qJ(4));
	t258 = qJD(1) * t240;
	t247 = qJD(5) + t258;
	t260 = t239 * t247;
	t241 = cos(qJ(5));
	t259 = t241 * t247;
	t242 = cos(qJ(4));
	t257 = qJD(1) * t242;
	t256 = qJD(4) * t240;
	t255 = qJD(4) * t242;
	t254 = qJD(5) * t240;
	t253 = qJD(5) * t242;
	t252 = t239 * t257;
	t251 = t241 * t257;
	t250 = t239 * t255;
	t249 = t241 * t255;
	t248 = -qJD(1) - t254;
	t246 = t239 * t253 + t241 * t256;
	t245 = t239 * t256 - t241 * t253;
	t244 = t248 * t241 - t250;
	t243 = t248 * t239 + t249;
	t238 = qJ(1) + pkin(8);
	t237 = cos(t238);
	t236 = sin(t238);
	t235 = t243 * t236 + t237 * t259;
	t234 = t244 * t236 - t237 * t260;
	t233 = -t236 * t259 + t243 * t237;
	t232 = t236 * t260 + t244 * t237;
	t1 = [t233, 0, 0, -t246 * t236 + t237 * t251, t234; t235, 0, 0, t236 * t251 + t246 * t237, -t232; 0, 0, 0, t239 * t254 - t249, t245; t232, 0, 0, t245 * t236 - t237 * t252, -t235; t234, 0, 0, -t236 * t252 - t245 * t237, t233; 0, 0, 0, t241 * t254 + t250, t246; t236 * t257 + t237 * t256, 0, 0, t236 * t255 + t237 * t258, 0; t236 * t256 - t237 * t257, 0, 0, t236 * t258 - t237 * t255, 0; 0, 0, 0, -t256, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end