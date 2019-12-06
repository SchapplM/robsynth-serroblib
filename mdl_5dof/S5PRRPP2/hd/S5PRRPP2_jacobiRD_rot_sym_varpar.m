% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:56
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(7));
	t15 = sin(pkin(7));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:56
	% EndTime: 2019-12-05 16:10:56
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
	t152 = sin(pkin(7));
	t161 = t152 * t165;
	t153 = cos(pkin(7));
	t160 = t153 * t165;
	t159 = t154 * t163 - t156 * t164;
	t158 = t154 * t164 + t156 * t163;
	t1 = [0, t159 * t153, t154 * t160 + (-t152 * t154 - t153 * t166) * qJD(3), 0, 0; 0, t159 * t152, t154 * t161 + (-t152 * t166 + t153 * t154) * qJD(3), 0, 0; 0, -t154 * t162 - t156 * t165, -t158, 0, 0; 0, t158 * t153, t156 * t160 + (-t152 * t156 + t153 * t167) * qJD(3), 0, 0; 0, t158 * t152, t156 * t161 + (t152 * t167 + t153 * t156) * qJD(3), 0, 0; 0, t154 * t165 - t156 * t162, t159, 0, 0; 0, -t160, 0, 0, 0; 0, -t161, 0, 0, 0; 0, t164, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:56
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t170 = sin(pkin(7));
	t173 = cos(qJ(2));
	t183 = t170 * t173;
	t171 = cos(pkin(7));
	t182 = t171 * t173;
	t172 = sin(qJ(2));
	t181 = qJD(2) * t172;
	t180 = qJD(2) * t173;
	t179 = qJD(3) * t172;
	t178 = qJD(3) * t173;
	t177 = t170 * t181;
	t176 = t171 * t181;
	t169 = qJ(3) + pkin(8);
	t167 = sin(t169);
	t168 = cos(t169);
	t175 = t167 * t179 - t168 * t180;
	t174 = t167 * t180 + t168 * t179;
	t1 = [0, t175 * t171, t167 * t176 + (-t167 * t170 - t168 * t182) * qJD(3), 0, 0; 0, t175 * t170, t167 * t177 + (t167 * t171 - t168 * t183) * qJD(3), 0, 0; 0, -t167 * t178 - t168 * t181, -t174, 0, 0; 0, t174 * t171, t168 * t176 + (t167 * t182 - t168 * t170) * qJD(3), 0, 0; 0, t174 * t170, t168 * t177 + (t167 * t183 + t168 * t171) * qJD(3), 0, 0; 0, t167 * t181 - t168 * t178, t175, 0, 0; 0, -t176, 0, 0, 0; 0, -t177, 0, 0, 0; 0, t180, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:57
	% EndTime: 2019-12-05 16:10:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t228 = sin(pkin(7));
	t231 = cos(qJ(2));
	t241 = t228 * t231;
	t229 = cos(pkin(7));
	t240 = t229 * t231;
	t230 = sin(qJ(2));
	t239 = qJD(2) * t230;
	t238 = qJD(2) * t231;
	t237 = qJD(3) * t230;
	t236 = qJD(3) * t231;
	t235 = t228 * t239;
	t234 = t229 * t239;
	t227 = qJ(3) + pkin(8);
	t225 = sin(t227);
	t226 = cos(t227);
	t233 = -t225 * t237 + t226 * t238;
	t232 = -t225 * t238 - t226 * t237;
	t1 = [0, -t233 * t229, t225 * t234 + (-t225 * t228 - t226 * t240) * qJD(3), 0, 0; 0, -t233 * t228, t225 * t235 + (t225 * t229 - t226 * t241) * qJD(3), 0, 0; 0, -t225 * t236 - t226 * t239, t232, 0, 0; 0, -t234, 0, 0, 0; 0, -t235, 0, 0, 0; 0, t238, 0, 0, 0; 0, t232 * t229, -t226 * t234 + (-t225 * t240 + t226 * t228) * qJD(3), 0, 0; 0, t232 * t228, -t226 * t235 + (-t225 * t241 - t226 * t229) * qJD(3), 0, 0; 0, -t225 * t239 + t226 * t236, t233, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end