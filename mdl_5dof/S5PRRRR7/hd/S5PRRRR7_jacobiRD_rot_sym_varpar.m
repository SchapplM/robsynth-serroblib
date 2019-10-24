% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:37
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:42
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(9));
	t15 = sin(pkin(9));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.05s
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
	t152 = sin(pkin(9));
	t161 = t152 * t165;
	t153 = cos(pkin(9));
	t160 = t153 * t165;
	t159 = t154 * t163 - t156 * t164;
	t158 = t154 * t164 + t156 * t163;
	t1 = [0, t159 * t153, t154 * t160 + (-t152 * t154 - t153 * t166) * qJD(3), 0, 0; 0, t159 * t152, t154 * t161 + (-t152 * t166 + t153 * t154) * qJD(3), 0, 0; 0, -t154 * t162 - t156 * t165, -t158, 0, 0; 0, t158 * t153, t156 * t160 + (-t152 * t156 + t153 * t167) * qJD(3), 0, 0; 0, t158 * t152, t156 * t161 + (t152 * t167 + t153 * t156) * qJD(3), 0, 0; 0, t154 * t165 - t156 * t162, t159, 0, 0; 0, -t160, 0, 0, 0; 0, -t161, 0, 0, 0; 0, t164, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (98->16), mult. (117->28), div. (0->0), fcn. (117->6), ass. (0->25)
	t211 = qJD(3) + qJD(4);
	t215 = sin(qJ(2));
	t227 = t211 * t215;
	t216 = cos(qJ(2));
	t226 = t211 * t216;
	t225 = qJD(2) * t215;
	t224 = qJD(2) * t216;
	t212 = qJ(3) + qJ(4);
	t209 = sin(t212);
	t223 = t209 * t226;
	t210 = cos(t212);
	t222 = t210 * t226;
	t213 = sin(pkin(9));
	t221 = t213 * t225;
	t214 = cos(pkin(9));
	t220 = t214 * t225;
	t219 = -t211 * t213 + t220;
	t218 = t211 * t214 + t221;
	t208 = t209 * t227 - t210 * t224;
	t217 = t209 * t224 + t210 * t227;
	t206 = t219 * t210 + t214 * t223;
	t205 = t219 * t209 - t214 * t222;
	t204 = t218 * t210 + t213 * t223;
	t203 = t218 * t209 - t213 * t222;
	t1 = [0, t208 * t214, t205, t205, 0; 0, t208 * t213, t203, t203, 0; 0, -t210 * t225 - t223, -t217, -t217, 0; 0, t217 * t214, t206, t206, 0; 0, t217 * t213, t204, t204, 0; 0, t209 * t225 - t222, t208, t208, 0; 0, -t220, 0, 0, 0; 0, -t221, 0, 0, 0; 0, t224, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (230->17), mult. (157->28), div. (0->0), fcn. (157->6), ass. (0->25)
	t229 = qJD(3) + qJD(4) + qJD(5);
	t233 = sin(qJ(2));
	t245 = t229 * t233;
	t234 = cos(qJ(2));
	t244 = t229 * t234;
	t243 = qJD(2) * t233;
	t242 = qJD(2) * t234;
	t230 = qJ(3) + qJ(4) + qJ(5);
	t227 = sin(t230);
	t241 = t227 * t244;
	t228 = cos(t230);
	t240 = t228 * t244;
	t231 = sin(pkin(9));
	t239 = t231 * t243;
	t232 = cos(pkin(9));
	t238 = t232 * t243;
	t237 = -t229 * t231 + t238;
	t236 = t229 * t232 + t239;
	t226 = t227 * t245 - t228 * t242;
	t235 = t227 * t242 + t228 * t245;
	t224 = t237 * t228 + t232 * t241;
	t223 = t237 * t227 - t232 * t240;
	t222 = t236 * t228 + t231 * t241;
	t221 = t236 * t227 - t231 * t240;
	t1 = [0, t226 * t232, t223, t223, t223; 0, t226 * t231, t221, t221, t221; 0, -t228 * t243 - t241, -t235, -t235, -t235; 0, t235 * t232, t224, t224, t224; 0, t235 * t231, t222, t222, t222; 0, t227 * t243 - t240, t226, t226, t226; 0, -t238, 0, 0, 0; 0, -t239, 0, 0, 0; 0, t242, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end