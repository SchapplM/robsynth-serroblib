% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR7_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t86 = sin(pkin(11));
t87 = sin(pkin(6));
t97 = t86 * t87;
t88 = cos(pkin(11));
t96 = t88 * t87;
t89 = cos(pkin(6));
t91 = sin(qJ(2));
t95 = t89 * t91;
t93 = cos(qJ(2));
t94 = t89 * t93;
t92 = cos(qJ(3));
t90 = sin(qJ(3));
t1 = [0, t97, t86 * t94 + t88 * t91, 0 (-t86 * t95 + t88 * t93) * t92 + t90 * t97, 0; 0, -t96, t86 * t91 - t88 * t94, 0 (t86 * t93 + t88 * t95) * t92 - t90 * t96, 0; 0, t89, -t87 * t93, 0, t87 * t91 * t92 + t89 * t90, 0;];
Jg_rot  = t1;
