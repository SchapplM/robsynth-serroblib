% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t23 = qJD(1) + qJD(2);
	t24 = qJ(1) + qJ(2);
	t25 = t23 * sin(t24);
	t20 = t23 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, 0, 0, 0; -t25, -t25, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->6), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t83 = qJD(1) + qJD(2);
	t90 = t83 * sin(pkin(9));
	t89 = t83 * cos(pkin(9));
	t84 = qJ(1) + qJ(2);
	t81 = sin(t84);
	t88 = t81 * t89;
	t82 = cos(t84);
	t87 = t82 * t90;
	t80 = t83 * t82;
	t79 = t83 * t81;
	t78 = t82 * t89;
	t77 = t81 * t90;
	t1 = [0, 0, 0, 0, 0; -t88, -t88, 0, 0, 0; t78, t78, 0, 0, 0; 0, 0, 0, 0, 0; t77, t77, 0, 0, 0; -t87, -t87, 0, 0, 0; 0, 0, 0, 0, 0; t80, t80, 0, 0, 0; t79, t79, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (119->15), mult. (132->24), div. (0->0), fcn. (132->6), ass. (0->23)
	t187 = qJD(1) + qJD(2);
	t189 = sin(pkin(9));
	t201 = t187 * t189;
	t190 = cos(pkin(9));
	t191 = sin(qJ(4));
	t200 = t190 * t191;
	t192 = cos(qJ(4));
	t199 = t190 * t192;
	t198 = qJD(4) * t189;
	t188 = qJ(1) + qJ(2);
	t185 = sin(t188);
	t197 = t185 * t201;
	t186 = cos(t188);
	t196 = t185 * t191 + t186 * t199;
	t195 = t185 * t192 - t186 * t200;
	t194 = -t185 * t199 + t186 * t191;
	t193 = t185 * t200 + t186 * t192;
	t184 = t186 * t201;
	t183 = -t193 * qJD(4) + t196 * t187;
	t182 = t194 * qJD(4) + t195 * t187;
	t181 = t195 * qJD(4) + t194 * t187;
	t180 = -t196 * qJD(4) + t193 * t187;
	t1 = [0, 0, 0, -t192 * t198, 0; t181, t181, 0, t182, 0; t183, t183, 0, -t180, 0; 0, 0, 0, t191 * t198, 0; t180, t180, 0, -t183, 0; t182, t182, 0, t181, 0; 0, 0, 0, 0, 0; -t197, -t197, 0, 0, 0; t184, t184, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (260->22), mult. (176->24), div. (0->0), fcn. (176->6), ass. (0->29)
	t239 = qJ(4) + qJ(5);
	t233 = sin(t239);
	t240 = qJ(1) + qJ(2);
	t234 = sin(t240);
	t254 = t233 * t234;
	t236 = cos(t240);
	t253 = t233 * t236;
	t235 = cos(t239);
	t252 = t234 * t235;
	t251 = t235 * t236;
	t237 = qJD(4) + qJD(5);
	t241 = sin(pkin(9));
	t250 = t237 * t241;
	t238 = qJD(1) + qJD(2);
	t249 = t238 * t241;
	t248 = t237 * t254;
	t247 = t237 * t251;
	t246 = t234 * t249;
	t245 = t235 * t250;
	t242 = cos(pkin(9));
	t244 = -t238 * t242 + t237;
	t243 = -t237 * t242 + t238;
	t232 = t236 * t249;
	t231 = t233 * t250;
	t226 = -t242 * t248 - t247 + (t242 * t251 + t254) * t238;
	t225 = t243 * t252 + t244 * t253;
	t224 = t243 * t253 + t244 * t252;
	t223 = -t242 * t247 - t248 + (t242 * t254 + t251) * t238;
	t1 = [0, 0, 0, -t245, -t245; t224, t224, 0, t225, t225; t226, t226, 0, -t223, -t223; 0, 0, 0, t231, t231; t223, t223, 0, -t226, -t226; t225, t225, 0, t224, t224; 0, 0, 0, 0, 0; -t246, -t246, 0, 0, 0; t232, t232, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end