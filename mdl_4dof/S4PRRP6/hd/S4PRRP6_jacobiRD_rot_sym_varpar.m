% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PRRP6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(6));
	t15 = sin(pkin(6));
	t1 = [0, -t16 * t19, 0, 0; 0, -t15 * t19, 0, 0; 0, -t20, 0, 0; 0, t16 * t20, 0, 0; 0, t15 * t20, 0, 0; 0, -t19, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:14
	% EndTime: 2019-12-29 12:22:14
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (18->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t154 = sin(qJ(3));
	t157 = cos(qJ(2));
	t167 = t154 * t157;
	t156 = cos(qJ(3));
	t166 = t156 * t157;
	t155 = sin(qJ(2));
	t165 = qJD(2) * t155;
	t164 = qJD(2) * t157;
	t163 = qJD(3) * t155;
	t162 = qJD(3) * t157;
	t152 = sin(pkin(6));
	t161 = t152 * t165;
	t153 = cos(pkin(6));
	t160 = t153 * t165;
	t159 = t154 * t163 - t156 * t164;
	t158 = t154 * t164 + t156 * t163;
	t1 = [0, t159 * t153, t154 * t160 + (-t152 * t154 - t153 * t166) * qJD(3), 0; 0, t159 * t152, t154 * t161 + (-t152 * t166 + t153 * t154) * qJD(3), 0; 0, -t154 * t162 - t156 * t165, -t158, 0; 0, t158 * t153, t156 * t160 + (-t152 * t156 + t153 * t167) * qJD(3), 0; 0, t158 * t152, t156 * t161 + (t152 * t167 + t153 * t156) * qJD(3), 0; 0, t154 * t165 - t156 * t162, t159, 0; 0, -t160, 0, 0; 0, -t161, 0, 0; 0, t164, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:15
	% EndTime: 2019-12-29 12:22:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (18->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t208 = sin(qJ(3));
	t211 = cos(qJ(2));
	t221 = t208 * t211;
	t210 = cos(qJ(3));
	t220 = t210 * t211;
	t209 = sin(qJ(2));
	t219 = qJD(2) * t209;
	t218 = qJD(2) * t211;
	t217 = qJD(3) * t209;
	t216 = qJD(3) * t211;
	t206 = sin(pkin(6));
	t215 = t206 * t219;
	t207 = cos(pkin(6));
	t214 = t207 * t219;
	t213 = -t208 * t217 + t210 * t218;
	t212 = -t208 * t218 - t210 * t217;
	t1 = [0, -t213 * t207, t208 * t214 + (-t206 * t208 - t207 * t220) * qJD(3), 0; 0, -t213 * t206, t208 * t215 + (-t206 * t220 + t207 * t208) * qJD(3), 0; 0, -t208 * t216 - t210 * t219, t212, 0; 0, -t214, 0, 0; 0, -t215, 0, 0; 0, t218, 0, 0; 0, t212 * t207, -t210 * t214 + (t206 * t210 - t207 * t221) * qJD(3), 0; 0, t212 * t206, -t210 * t215 + (-t206 * t221 - t207 * t210) * qJD(3), 0; 0, -t208 * t219 + t210 * t216, t213, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end