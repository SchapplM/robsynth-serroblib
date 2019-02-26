% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t84 = sin(pkin(10));
t85 = sin(pkin(6));
t95 = t84 * t85;
t86 = cos(pkin(10));
t94 = t86 * t85;
t87 = cos(pkin(6));
t89 = sin(qJ(2));
t93 = t87 * t89;
t91 = cos(qJ(2));
t92 = t87 * t91;
t90 = cos(qJ(4));
t88 = sin(qJ(4));
t1 = [0, t95, 0, -t84 * t93 + t86 * t91, 0, t90 * t95 + (t84 * t92 + t86 * t89) * t88; 0, -t94, 0, t84 * t91 + t86 * t93, 0, -t90 * t94 + (t84 * t89 - t86 * t92) * t88; 0, t87, 0, t85 * t89, 0, -t85 * t91 * t88 + t87 * t90;];
Jg_rot  = t1;
