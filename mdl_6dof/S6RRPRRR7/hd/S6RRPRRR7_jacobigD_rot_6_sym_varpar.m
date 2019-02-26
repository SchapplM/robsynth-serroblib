% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:35
% EndTime: 2019-02-26 21:57:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (28->7), mult. (92->12), div. (0->0), fcn. (92->6), ass. (0->15)
t171 = qJD(2) - qJD(4);
t162 = sin(qJ(1));
t170 = qJD(1) * t162;
t165 = cos(qJ(1));
t169 = qJD(1) * t165;
t160 = sin(qJ(4));
t161 = sin(qJ(2));
t163 = cos(qJ(4));
t164 = cos(qJ(2));
t168 = t160 * t164 - t161 * t163;
t166 = t171 * (t160 * t161 + t163 * t164);
t159 = t171 * t168;
t158 = -t166 * t162 + t168 * t169;
t157 = -t166 * t165 - t168 * t170;
t1 = [0, t169, 0, -t169, t157, t157; 0, t170, 0, -t170, t158, t158; 0, 0, 0, 0, t159, t159;];
JgD_rot  = t1;
