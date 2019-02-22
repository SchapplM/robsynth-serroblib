% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:32:26
% EndTime: 2019-02-22 10:32:26
% DurationCPUTime: 0.04s
% Computational Cost: add. (83->14), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t91 = qJ(5) + qJ(6);
t88 = sin(t91);
t92 = sin(qJ(3));
t96 = t92 * t88;
t89 = cos(t91);
t95 = t92 * t89;
t93 = cos(qJ(3));
t85 = t93 * t88;
t94 = t93 * t89;
t90 = qJ(1) + pkin(10);
t87 = cos(t90);
t86 = sin(t90);
t84 = -t86 * t96 + t87 * t89;
t83 = t86 * t95 + t87 * t88;
t82 = t86 * t89 + t87 * t96;
t81 = -t86 * t88 + t87 * t95;
t1 = [t84, 0, t87 * t85, 0, t81, t81; t82, 0, t86 * t85, 0, t83, t83; 0, 0, t96, 0, -t94, -t94; -t83, 0, t87 * t94, 0, -t82, -t82; t81, 0, t86 * t94, 0, t84, t84; 0, 0, t95, 0, t85, t85; -t86 * t93, 0, -t87 * t92, 0, 0, 0; t87 * t93, 0, -t86 * t92, 0, 0, 0; 0, 0, t93, 0, 0, 0;];
JR_rot  = t1;
