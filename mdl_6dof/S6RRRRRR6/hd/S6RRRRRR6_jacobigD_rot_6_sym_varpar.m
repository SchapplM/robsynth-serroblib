% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:05
% EndTime: 2019-02-26 22:50:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (78->18), mult. (152->40), div. (0->0), fcn. (156->8), ass. (0->30)
t220 = sin(pkin(6));
t223 = sin(qJ(1));
t238 = t220 * t223;
t225 = cos(qJ(1));
t237 = t220 * t225;
t222 = sin(qJ(2));
t236 = t222 * t223;
t235 = t222 * t225;
t224 = cos(qJ(2));
t234 = t223 * t224;
t233 = t225 * t224;
t232 = qJD(1) * t220;
t219 = qJ(3) + qJ(4);
t216 = sin(t219);
t231 = qJD(2) * t216;
t230 = qJD(2) * t220;
t221 = cos(pkin(6));
t229 = t221 * t233 - t236;
t228 = t221 * t234 + t235;
t227 = t221 * t235 + t234;
t226 = -t221 * t236 + t233;
t218 = qJD(3) + qJD(4);
t217 = cos(t219);
t215 = t222 * t230;
t214 = t228 * qJD(1) + t227 * qJD(2);
t213 = t229 * qJD(1) + t226 * qJD(2);
t212 = t220 * t222 * t218 * t217 + (t218 * t221 + t224 * t230) * t216;
t211 = (-t216 * t237 + t227 * t217) * t218 + t229 * t231 + (t226 * t216 - t217 * t238) * qJD(1);
t210 = (t216 * t238 + t226 * t217) * t218 - t228 * t231 + (-t227 * t216 - t217 * t237) * qJD(1);
t1 = [0, t225 * t232, t213, t213, t210, t210; 0, t223 * t232, t214, t214, t211, t211; 0, 0, t215, t215, t212, t212;];
JgD_rot  = t1;
