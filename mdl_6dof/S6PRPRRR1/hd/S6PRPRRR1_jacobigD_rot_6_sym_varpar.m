% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:44
% EndTime: 2019-02-26 19:53:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (49->15), mult. (127->33), div. (0->0), fcn. (134->10), ass. (0->23)
t216 = sin(pkin(12));
t219 = cos(pkin(12));
t222 = sin(qJ(2));
t223 = cos(qJ(2));
t225 = t222 * t216 - t223 * t219;
t229 = t225 * qJD(2);
t226 = t216 * t223 + t219 * t222;
t210 = t226 * qJD(2);
t214 = qJD(4) + qJD(5);
t215 = qJ(4) + qJ(5);
t228 = cos(t215) * t214;
t218 = sin(pkin(6));
t221 = cos(pkin(6));
t227 = t214 * t218 + t221 * t229;
t220 = cos(pkin(11));
t217 = sin(pkin(11));
t212 = sin(t215);
t208 = t226 * t221;
t206 = t221 * t210;
t205 = t218 * t210;
t204 = -t217 * t206 - t220 * t229;
t203 = t220 * t206 - t217 * t229;
t1 = [0, 0, 0, t204, t204 (-t217 * t208 - t220 * t225) * t228 + (-t220 * t210 + t217 * t227) * t212; 0, 0, 0, t203, t203 (t220 * t208 - t217 * t225) * t228 + (-t217 * t210 - t227 * t220) * t212; 0, 0, 0, t205, t205, t221 * t214 * t212 + (-t212 * t229 + t226 * t228) * t218;];
JgD_rot  = t1;
