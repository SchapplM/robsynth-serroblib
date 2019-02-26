% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR13_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobigD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
t192 = sin(pkin(6));
t205 = t192 * cos(pkin(7));
t195 = sin(qJ(2));
t196 = sin(qJ(1));
t204 = t195 * t196;
t198 = cos(qJ(1));
t203 = t195 * t198;
t197 = cos(qJ(2));
t202 = t196 * t197;
t201 = t197 * t198;
t200 = qJD(1) * t192;
t191 = sin(pkin(7));
t199 = qJD(2) * t191;
t194 = cos(pkin(6));
t1 = [0, t198 * t200 -(t194 * t204 - t201) * t199 + (-(-t194 * t201 + t204) * t191 + t198 * t205) * qJD(1), 0, 0, 0; 0, t196 * t200 -(-t194 * t203 - t202) * t199 + (-(-t194 * t202 - t203) * t191 + t196 * t205) * qJD(1), 0, 0, 0; 0, 0, t192 * t195 * t199, 0, 0, 0;];
JgD_rot  = t1;
