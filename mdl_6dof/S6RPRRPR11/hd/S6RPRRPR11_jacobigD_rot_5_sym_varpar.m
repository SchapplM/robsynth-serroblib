% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:40
% EndTime: 2019-02-26 21:06:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
t225 = sin(pkin(7));
t226 = sin(pkin(6));
t246 = t226 * t225;
t228 = cos(pkin(7));
t232 = cos(qJ(3));
t245 = t228 * t232;
t224 = sin(pkin(12));
t231 = sin(qJ(1));
t244 = t231 * t224;
t227 = cos(pkin(12));
t243 = t231 * t227;
t233 = cos(qJ(1));
t242 = t233 * t224;
t241 = t233 * t227;
t240 = t231 * t246;
t239 = t233 * t246;
t238 = qJD(1) * t226 * t228;
t229 = cos(pkin(6));
t237 = t229 * t241 - t244;
t236 = -t229 * t243 - t242;
t235 = t229 * t242 + t243;
t234 = -t229 * t244 + t241;
t230 = sin(qJ(3));
t223 = t236 * qJD(1);
t222 = t237 * qJD(1);
t1 = [0, 0, t222 * t225 + t233 * t238, t222 * t245 + (t234 * t232 + (t236 * t228 + t240) * t230) * qJD(3) + (-t235 * t230 - t232 * t239) * qJD(1), 0, 0; 0, 0, -t223 * t225 + t231 * t238, -t223 * t245 + (t235 * t232 + (t237 * t228 - t239) * t230) * qJD(3) + (t234 * t230 - t232 * t240) * qJD(1), 0, 0; 0, 0, 0 (t225 * t229 * t230 + (t227 * t228 * t230 + t224 * t232) * t226) * qJD(3), 0, 0;];
JgD_rot  = t1;
