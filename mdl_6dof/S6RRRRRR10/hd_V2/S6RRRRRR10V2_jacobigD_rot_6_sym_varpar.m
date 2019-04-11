% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10V2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (77->28), mult. (112->57), div. (0->0), fcn. (114->8), ass. (0->26)
t232 = qJ(2) + qJ(3);
t230 = cos(t232);
t237 = cos(qJ(4));
t251 = t230 * t237;
t231 = qJD(2) + qJD(3);
t250 = t231 * t230;
t235 = sin(qJ(1));
t249 = t231 * t235;
t248 = t231 * t237;
t233 = sin(qJ(5));
t247 = t233 * t237;
t238 = cos(qJ(1));
t246 = t238 * t231;
t227 = qJD(1) * t235;
t228 = qJD(1) * t238;
t234 = sin(qJ(4));
t245 = qJD(4) * t234;
t244 = -qJD(5) + t248;
t243 = qJD(4) * t230 - qJD(1);
t242 = qJD(1) * t230 - qJD(4);
t241 = t243 * t237;
t240 = t242 * t235;
t229 = sin(t232);
t239 = -t229 * t227 + t230 * t246;
t236 = cos(qJ(5));
t1 = [0, t228, t228, t239, t238 * t241 + (-t229 * t246 - t240) * t234, -t240 * t247 + (-t244 * t229 - t243 * t234) * t233 * t238 + ((t235 * t234 + t238 * t251) * qJD(5) - t239) * t236; 0, t227, t227, t229 * t228 + t230 * t249, t235 * t241 + (-t229 * t249 + t242 * t238) * t234 (t242 * t247 + (-qJD(1) * t229 - t234 * qJD(5)) * t236) * t238 + ((qJD(1) * t234 - t229 * t248 - t230 * t245) * t233 - t236 * t250 + (t229 * t233 + t236 * t251) * qJD(5)) * t235; 0, 0, 0, t231 * t229, qJD(4) * t229 * t237 + t234 * t250, t244 * t233 * t230 + (-t233 * t245 + (qJD(5) * t237 - t231) * t236) * t229;];
JgD_rot  = t1;
