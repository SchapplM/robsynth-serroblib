% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:45
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
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
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t16 = qJD(1) * sin(qJ(1));
	t15 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; t15, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t15, 0, 0, 0, 0, 0; -t16, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:45
	% EndTime: 2019-10-09 23:42:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(4));
	t40 = qJD(4) * t34;
	t36 = cos(qJ(4));
	t39 = qJD(4) * t36;
	t38 = qJD(4) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = -t34 * t41 - t35 * t39;
	t31 = -t34 * t38 - t36 * t42;
	t30 = t34 * t42 - t36 * t38;
	t1 = [t32, 0, 0, t31, 0, 0; -t30, 0, 0, t33, 0, 0; 0, 0, 0, -t39, 0, 0; -t33, 0, 0, t30, 0, 0; t31, 0, 0, t32, 0, 0; 0, 0, 0, t40, 0, 0; t42, 0, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:45
	% EndTime: 2019-10-09 23:42:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t160 = sin(qJ(1));
	t167 = qJD(1) * t160;
	t162 = cos(qJ(1));
	t166 = qJD(1) * t162;
	t159 = sin(qJ(4));
	t165 = qJD(4) * t159;
	t161 = cos(qJ(4));
	t164 = qJD(4) * t161;
	t163 = qJD(4) * t162;
	t158 = -t160 * t165 + t161 * t166;
	t157 = t159 * t166 + t160 * t164;
	t156 = t159 * t163 + t161 * t167;
	t155 = t159 * t167 - t161 * t163;
	t1 = [t167, 0, 0, 0, 0, 0; -t166, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t157, 0, 0, t156, 0, 0; t155, 0, 0, -t158, 0, 0; 0, 0, 0, t164, 0, 0; t158, 0, 0, -t155, 0, 0; t156, 0, 0, t157, 0, 0; 0, 0, 0, -t165, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:46
	% EndTime: 2019-10-09 23:42:46
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t236 = cos(qJ(1));
	t235 = cos(qJ(4));
	t255 = qJD(1) * t235;
	t240 = qJD(6) + t255;
	t258 = t240 * t236;
	t233 = sin(qJ(1));
	t232 = sin(qJ(4));
	t251 = qJD(4) * t236;
	t243 = t232 * t251;
	t257 = t240 * t233 + t243;
	t256 = qJD(1) * t233;
	t254 = qJD(1) * t236;
	t253 = qJD(4) * t232;
	t252 = qJD(4) * t235;
	t231 = sin(qJ(6));
	t250 = qJD(6) * t231;
	t249 = qJD(6) * t232;
	t248 = qJD(6) * t235;
	t234 = cos(qJ(6));
	t247 = t234 * t253;
	t246 = t234 * t249;
	t245 = t233 * t253;
	t244 = t233 * t252;
	t242 = t235 * t251;
	t241 = qJD(1) + t248;
	t239 = t241 * t236;
	t238 = t232 * t254 + t244;
	t237 = -t232 * t256 + t242;
	t230 = -t257 * t231 + t234 * t239;
	t229 = t231 * t239 + t257 * t234;
	t228 = t241 * t234 * t233 + (-t245 + t258) * t231;
	t227 = -t234 * t258 + (t241 * t231 + t247) * t233;
	t1 = [t228, 0, 0, t237 * t231 + t236 * t246, 0, t229; -t230, 0, 0, t238 * t231 + t233 * t246, 0, t227; 0, 0, 0, -t231 * t253 + t234 * t248, 0, -t231 * t249 + t234 * t252; -t227, 0, 0, t234 * t242 + (-t234 * t256 - t236 * t250) * t232, 0, t230; t229, 0, 0, t234 * t244 + (-t233 * t250 + t234 * t254) * t232, 0, t228; 0, 0, 0, -t231 * t248 - t247, 0, -t231 * t252 - t246; -t238, 0, 0, -t233 * t255 - t243, 0, 0; t237, 0, 0, t235 * t254 - t245, 0, 0; 0, 0, 0, -t252, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end