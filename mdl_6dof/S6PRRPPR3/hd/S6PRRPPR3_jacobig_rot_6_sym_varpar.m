% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPPR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:24
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t85 = sin(pkin(10));
t86 = sin(pkin(6));
t96 = t85 * t86;
t87 = cos(pkin(10));
t95 = t87 * t86;
t88 = cos(pkin(6));
t90 = sin(qJ(2));
t94 = t88 * t90;
t92 = cos(qJ(2));
t93 = t88 * t92;
t91 = cos(qJ(3));
t89 = sin(qJ(3));
t1 = [0, t96, t85 * t93 + t87 * t90, 0, 0 (-t85 * t94 + t87 * t92) * t91 + t89 * t96; 0, -t95, t85 * t90 - t87 * t93, 0, 0 (t85 * t92 + t87 * t94) * t91 - t89 * t95; 0, t88, -t86 * t92, 0, 0, t86 * t90 * t91 + t88 * t89;];
Jg_rot  = t1;
