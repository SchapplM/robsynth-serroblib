% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRPRP4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->5), mult. (10->8), div. (0->0), fcn. (10->4), ass. (0->6)
	t26 = qJD(2) * sin(pkin(7));
	t25 = qJD(2) * cos(pkin(7));
	t22 = qJ(2) + pkin(8);
	t21 = cos(t22);
	t20 = sin(t22);
	t1 = [0, -t21 * t25, 0, 0, 0; 0, -t21 * t26, 0, 0, 0; 0, -qJD(2) * t20, 0, 0, 0; 0, t20 * t25, 0, 0, 0; 0, t20 * t26, 0, 0, 0; 0, -qJD(2) * t21, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->16), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t168 = sin(pkin(7));
	t170 = sin(qJ(4));
	t184 = t168 * t170;
	t171 = cos(qJ(4));
	t183 = t168 * t171;
	t169 = cos(pkin(7));
	t182 = t169 * t170;
	t181 = t169 * t171;
	t167 = qJ(2) + pkin(8);
	t165 = sin(t167);
	t180 = qJD(2) * t165;
	t179 = qJD(2) * t170;
	t178 = qJD(2) * t171;
	t177 = qJD(4) * t170;
	t176 = qJD(4) * t171;
	t175 = t168 * t180;
	t174 = t169 * t180;
	t166 = cos(t167);
	t173 = t165 * t177 - t166 * t178;
	t172 = t165 * t176 + t166 * t179;
	t1 = [0, t173 * t169, 0, t170 * t174 + (-t166 * t181 - t184) * qJD(4), 0; 0, t173 * t168, 0, t170 * t175 + (-t166 * t183 + t182) * qJD(4), 0; 0, -t165 * t178 - t166 * t177, 0, -t172, 0; 0, t172 * t169, 0, t171 * t174 + (t166 * t182 - t183) * qJD(4), 0; 0, t172 * t168, 0, t171 * t175 + (t166 * t184 + t181) * qJD(4), 0; 0, t165 * t179 - t166 * t176, 0, t173, 0; 0, -t174, 0, 0, 0; 0, -t175, 0, 0, 0; 0, qJD(2) * t166, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:32
	% EndTime: 2019-12-05 15:36:32
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (45->17), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->21)
	t227 = sin(pkin(7));
	t229 = sin(qJ(4));
	t243 = t227 * t229;
	t230 = cos(qJ(4));
	t242 = t227 * t230;
	t228 = cos(pkin(7));
	t241 = t228 * t229;
	t240 = t228 * t230;
	t226 = qJ(2) + pkin(8);
	t224 = sin(t226);
	t239 = qJD(2) * t224;
	t238 = qJD(2) * t229;
	t237 = qJD(2) * t230;
	t236 = qJD(4) * t229;
	t235 = qJD(4) * t230;
	t234 = t227 * t239;
	t233 = t228 * t239;
	t225 = cos(t226);
	t232 = -t224 * t236 + t225 * t237;
	t231 = -t224 * t235 - t225 * t238;
	t1 = [0, -t232 * t228, 0, t229 * t233 + (-t225 * t240 - t243) * qJD(4), 0; 0, -t232 * t227, 0, t229 * t234 + (-t225 * t242 + t241) * qJD(4), 0; 0, -t224 * t237 - t225 * t236, 0, t231, 0; 0, -t233, 0, 0, 0; 0, -t234, 0, 0, 0; 0, qJD(2) * t225, 0, 0, 0; 0, t231 * t228, 0, -t230 * t233 + (-t225 * t241 + t242) * qJD(4), 0; 0, t231 * t227, 0, -t230 * t234 + (-t225 * t243 - t240) * qJD(4), 0; 0, -t224 * t238 + t225 * t235, 0, t232, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end